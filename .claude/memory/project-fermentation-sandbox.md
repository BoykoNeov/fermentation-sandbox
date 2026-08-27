---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-08-17T19:56:05.039Z
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation engine in Python (uv, scipy/numpy/pydantic). Repo: https://github.com/BoykoNeov/fermentation-sandbox (branch `main`).

**Session-boot context: PROHIBITIONS and POINTERS only** — not a changelog. Every bullet is *what it forbids* +
the D-record to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue past
it from this file.** **Detail lives in `prohibitions/` (per subject) and `lessons/` (epistemics), reached BY PATH —
NO `MEMORY.md` row, so they cost nothing until read (D-185; lessons split the same way 2026-08-26).**
**Caps: 8 lines per BLOCK here, 14 in `CLAUDE.md`, 320 BYTES per index row here and in `lessons/` (chars until
2026-08-26 — the loader counts bytes) — and NO whole-file total, removed at D-177. `MEMORY.md`'s headroom against
the harness load limit is REPORTED every write; it truncates SILENTLY, which is how it ran over twice.**
(`.claude/hooks/check_memory_size.py`; [[feedback-batch-end-ritual]]). Distil the NEW block; **never evict an old prohibition to buy a line.**

## Where the records are
- `docs/DECISIONS.md` — canonical archive, **tens of thousands of lines: never read it linearly.** Generated top block gives a
  subsystem cut, ordered list, and **correction map (⚠)**. **The ⚠ lives only in the index** — check a record's
  index row before trusting it. Append per `CLAUDE.md`'s `Corrects:`/`Flags:`, then `tools/gen_decisions_toc.py`
  (**edit `TOPIC_RULES` if a new record buckets nowhere**). **File is LF.**
- `docs/ARCHITECTURE.md` — **the structure map: where things live, what imports what. Rebuilt from measurement
  at D-184** (it had gone 127 archive-commits stale). Counts are DERIVED — a snippet at its foot re-derives every
  one; if they disagree the code is right. **Update it in the SAME commit** as any Process/slot/param-file/package
  change. `docs/plans/*.md` are **FROZEN LOGS, bannered as such — never read them for what is built/open/next**;
  the "keep the plans updated" rule is **RETIRED** (D-184). `CLAUDE.md` = prime directives + archive conventions.

## Status (2026-08-27)
M0/M1/M2 **complete**. **M3** (sensory/OAV + Tier-3 aging, D-66) at **D-241** (D-231 is tooling, not model); sensory 1a/1b + **D-139's leftovers ALL closed** (D-148/9). Suite **1979 + 5 xfail** (3 D-188 Herzan; 2 D-215's own) — MEASURED off a green full run, never inferred: the ledger read 1861 for two beats.
**§2.2's beer criterion PASSES since D-223** (6.04 d in 5-7) but is partly self-referential and is labelled so in its own docstring. **"Blocked on external sourcing" wrong 6× (D-191/196/199/208/209/211)** — and **D-230 adds a 7th shape: the source was in THIS REPO**, transcribed at D-209 for another derivation.
**D-203/205/206 REFUSED** the sotolon-ascorbate route, Pham's pH+ethanol terms and the Strecker split — "expressible" ≠ "identifiable"; D-202 completed Fig 24.12's top group and D-204 shipped its pin.
Slot/Process/oxidative counts live in `docs/ARCHITECTURE.md` — never restate here, that rotted it (D-184).
Beer acid-base = **NINE** beats BUILT (D-178→D-183, D-207→D-209, D-211, **D-239**), `ACID_STATE` NOT medium-agnostic (D-179); **D-212 BUILT NOTHING** — day-1 pH admits an acetic WINDOW that Tyrell's own 145 sits OUTSIDE at all 3 arms. **D-239 BUILT D-209 §8's buffer-removal half and its day-7 cost puts the high `z̄` edge OUTSIDE Tyrell's envelope ON PURPOSE — never xfail that edge.**
**D-232 REFUSED BOTH ways of deciding D-230's residue** (a settling profile that flips SIGN across `mu_max`'s band; a pH clock that only re-measures a known +0.162 day-1 miss) and **ADDED a third branch** — it widened the ambiguity, never narrowed it.
**D-213 BUILT beer's wort O₂ — sourced and INERT**, and its decline of the O₂→growth coupling rested on "none of the six predictions is reachable in the default set": **growth EXTENT now is** (D-230 §7). Owner's call to re-open, never a beat's.
**D-234 RAN D-233's census (32 names) and found the repo in a state D-186's own docstring FORBADE**; the owner called BOTH repairs and they SHIPPED — **D-235** (mutations get the RUNNING param map, so `set_ph` re-anchors per member: 0.07896 → 1.98e-11; it also exposed `add_dap` reading z̄ off the compile map, 0.0108 pH, repaired too) and **D-236** (`y0_for_member` re-seeds `copper` per member: 16.65 % browning spread → bit-identical, CONDITIONAL on the scenario not naming `copper_gpl`). **The census has NO LIVE row left.**

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

## Live prohibitions — LEDGER (detail is one hop away, by path)
**A subject listed here is SETTLED: do not re-propose any of it as unbuilt.** That guarantee is
unconditional and lives here; the linked file only tells you *why*. **Read the file before proposing
work that touches its subject.** Split out at D-185 from 320 inline lines
[[feedback-a-doc-rots-where-it-duplicates]].
- **Sampling surfaces (D-153 → D-162)** — BOTH archive-wide sweeps DONE; four surfaces and two
  distributions PINNED; D-157's contradiction CLOSED (never re-narrow 0.084); `reads` has TWO masters;
  closure ordering SCOPED. Do not re-run, re-audit or "simplify" any of it.
  → `.claude/memory/prohibitions/sampling-surfaces.md`
- **Band edges & provenance (D-163 → D-176)** — the 73-arm edge sweep DONE; class (d), BOTH wide-band
  mechanisms, the switch-site census and the "buildable 55" all CLOSED; the O2 partition
  reparameterised; oak's share is a JOINT. **ZERO edges moved** bar `E_a_oak_extraction`. The Q10
  sweep is **NOT idempotent under its own repair — never re-run it expecting a defect count.**
  → `.claude/memory/prohibitions/band-edges-and-provenance.md`
- **Oxidation — direct set + burst (D-132 → D-137, D-147, D-149 → D-152)** — burst is **WIRED and
  NON-DEFAULT**, not dead code and not refused; copper 600 L/g CLOSED and BOUNDED; a pH term on
  activation is REFUSED; the O2 gate is Fe(II)+O2 with SO2 INVERTED in role.
  → `.claude/memory/prohibitions/oxidation-direct-and-burst.md`
- **Oxidation — cascade + SO2 dosing (D-141 → D-148)** — cascade is **BUILT and NON-DEFAULT**, benchmark
  ACTIVE not open work; must-dosing REFUSED; "the sim under-binds SO2" WITHDRAWN; §2.4 CLOSED.
  → `.claude/memory/prohibitions/oxidation-cascade-and-dosing.md`
- **Aroma compounds + Milestone-3 tail (D-117 → D-136, D-146, D-176)** — methionine sink is **BLOCKED
  not deferred**; **do NOT rebuild bottle reduction as thioacetate/disulfide precursors** (D-101's
  mechanism guess was WRONG) and sulfide release is **temperature-FLAT**; technical cork ships BELOW
  screwcap; **never put 0.963 in a sampled field**; the `oav`→`magnitude` rename must NOT be
  "finished"; EtOAc eq is Berthelot-coupled. → `.claude/memory/prohibitions/aroma-and-milestone-3-tail.md`
- **Aging-pH anchor (D-186)** — `set_ph` is **BUILT**: D-150's "no way to set an aging pH" is FALSE,
  though ~8 records' copy-forward lists still say it. **Cation-MOVING, never a pH dial**; `add_acid`
  is unchanged; the reachability check **cannot** move to compile; the opt-in gate's reason DIFFERS
  per medium; adds **no** pH-rate dependence. → `.claude/memory/prohibitions/aging-ph-anchor.md`
- **Closure oxygen — steady + bottling burst (D-136 → D-162, D-187)** — BOTH columns of Lopes'
  table now ship: `seal_bottle` is **BUILT** and D-136's "not for lack of data" is spent. Dose is
  the first month **NET of steady** and must not precede `begin_aging`; screwcap is a **BOUND**,
  never a midpoint; the burst ordering is NOT the closure ordering; line oxygen stays with
  `add_oxygen`. OTR(T) and bottle format still blocked.
  → `.claude/memory/prohibitions/closure-oxygen.md`
- **Acetaldehyde ladder — Herzan 2020 (D-188)** — MEASURED, attributed, **NOT built** (only anchor is
  the table under test). **D-108's "no maturation source" is FALSE**; the `0-vs-2.7` floor is CLOSED as
  the closure menu's span. **Never a per-cell ratio off this ladder** (the model emits one constant);
  never group must-sulfited vs not; **D-47 protection WORKS** — competition owns the SIGN, a missing
  post-dryness source owns the SIZE. → `.claude/memory/prohibitions/acetaldehyde-ladder.md`
- **Keto-acid excretion SHAPE — option B (D-49 → D-195)** — MEASURED and **REFUSED for BOTH
  α-ketobutyrate (D-189) and PYRUVATE (D-195)**: a growth-linked source with the flux-linked sink
  **drains the pool it feeds**, arithmetic **general to all three pools**. **Never argue pyruvate's
  as invisible-gain — option B WORKS there.** Refused because draining the residual leaves
  **Miao's measured band**, and a control pins that on the **RESIDUAL, not the shape**. **No guard
  owed either time.** D-107's shape-diagnosis is **corrected**: the miss is a competition.
  → `.claude/memory/prohibitions/keto-acid-excretion-shape.md`
- **Carbonyl release + the binding constants (D-190)** — release is **ALREADY EMERGENT** for
  acetaldehyde (stateless equilibrium, never build it); the aroma carbonyls bind **NOTHING** so
  Bueno's ordering is inexpressible and methional is **BLOCKED on a Kd**; diacetyl REFUSED on the
  table's reliability; the pyruvate/α-KG ordering is corrected **in prose only — no value moved**,
  and its reversal is reachable in **2.0 %** of draws. → `.claude/memory/prohibitions/carbonyl-release-and-binding.md`
- **Residual copper after fining (D-191)** — the credit is **BUILT**: `add_copper` writes the
  `copper` slot, D-44's "CuS drops out with the lees" is **FALSE**, D-149's "two coppers never
  meet" **CLOSED**, and the docstring's *"nothing sources one"* was false when written. Never
  midpoint the 0.95; the **sulfide half is MEASURED, not built** (NOT a live candidate — spent);
  D-45's mercaptide carbon flow is **FLAGGED not fixed** — **both spent at D-193, see the row above**.
  → `.claude/memory/prohibitions/residual-copper-fining.md`
- **Fined sulfur's destination (D-193)** — the transfer into `bound_h2s`/`bound_methanethiol` is
  **BUILT**: fining complexes the sulfur, never destroys it, so a fining is **no longer permanent**.
  **1:1, WHOLE mass — the sulfur is NOT scaled by the 0.95** (that variant was built and REJECTED;
  copper 0.95 and sulfur 1.0 differ on purpose and a test forbids "harmonising" them). D-45's carbon
  booking is CORRECTED to **zero**, not shrunk. **D-135's rate-coupling refusal and its unmodelled
  SPONTANEOUS formation route both STAND**; its "a reservoir fining never touched" docstring does not.
  → `.claude/memory/prohibitions/fined-sulfur-destination.md`
- **Sweet-wine scenario anchor (D-194)** — `_SWEET_BRIX` is **38.0**, D-192's OPEN item is
  **CLOSED**. At 70 it was a must **still fermenting** (0.50 % ABV at the breakpoint), not a
  weak wine. **Never anchor on the residual** — 45 °Brix was REJECTED as hiding the 20.2 % ABV
  ceiling, which stays visible. `assert S > 50` passed for the WRONG REASON, and the green
  mutation is what OWED the guard. → `.claude/memory/prohibitions/sweet-wine-anchor.md`
- **Osmotic inhibition at high sugar (D-192)** — **BUILT**, wine-only, threshold **300 g/L =
  the TOP of Coleman's envelope, never the Handbook's printed 200** (that is inside the
  keystone's own fit). The Handbook's "300 yields less than 200" is **FORFEITED on measurement**;
  the Haldane form is **refuted on shape** (spans 1.51×, needs ~19×); never a hard zero (absorbing
  state); `K`/`n` are a **derived pair, not two bands**. `_SWEET_BRIX` is CLOSED at 38.0 (D-194).
  → `.claude/memory/prohibitions/osmotic-high-sugar.md`
- **Post-Fenton second O₂ (D-196)** — **BUILT** (cascade-only, still non-default): the radical
  takes an O₂, `share × activation`. **NEVER "unsourced"** — 20 records said so, the book was on
  disk. Per-ACETALDEHYDE not per-H₂O₂; **its "upper bound on the net draw" gloss is FALSE — D-198
  measured a LOWER bound**; Gate 1 not breached. Double-count **died on supply-limitation**, D-141's
  survives. The state-dependence is **load-bearing** — a flat surcharge of equal cost breaks BOTH
  asymptotes. Copper re-fit is **exactly** a no-op. D-142 **CONFIRMED, not corrected**.
  → `.claude/memory/prohibitions/post-fenton-second-o2.md`
- **Quinone double charge in the DEFAULT set (D-197)** — MEASURED: D-75's lump is **FIVE draws
  taking 36-76 % of the O₂**, not the ~2 % its own defence cites (that number is TRUE and about the
  **smallest** member; anthocyanin fading takes ~1000× more). **Never de-duplicate in the parallel
  frame** — LOAD-BEARING: 0.384 → 0.109 vs Carrascón's measured 0.88-1.25. The two-stage rework
  **is the cascade**, not unbuilt. **Use the 90 % limb, not the 7-day cap** (the cap flatters it to
  "15 %"; real gap ~2.3×, an accepted deviation). NO guard was owed; margins ride on free SO₂.
  → `.claude/memory/prohibitions/quinone-double-charge.md`
- **Hydroperoxyl recycling limb (D-198)** — **REFUSED, and now MEASURED not inherited.** D-197
  expected a "documented no"; that passage is the **initiation** node's HO₂· (p. 327-8), a
  DIFFERENT reaction from the ethanol-limb co-product (p. 332). Four grounds: no source for this
  limb's co-product; **three sourced fates disagree**; a **pole at s=1** (unsulfited) with no
  published quench ratio; and nothing here can **separate** the candidates (O₂ budget spans
  0.078 %). It is a **RATE-LAW** change (`A/(1−s)`), not stoichiometry. Do not re-open by
  re-sweeping the 24 texts. → `.claude/memory/prohibitions/hydroperoxyl-recycling.md`
- **Quinone branching (D-145 → D-202)** — the named pull **is ON DISK** and **still cannot close
  it** (ranks NUCLEOPHILES; the 84.78 % consumer has none). **Blocked on STRUCTURE**, but D-202
  **RULES OUT D-141's polymerisation-band closure**. **GSH is MEASURED, NOT BUILT** (0.32 %; never
  price it off Ch. 24's bound, ~10× loose vs Ch. 5). **1.7481 at ×0.01 REFUSED.** High edge 0.9738,
  **not D-174's**. Sulfonate is the **MINOR** product. **H2S BUILT (D-201), ascorbate BUILT
  (D-202) — the top group is COMPLETE**; D-200's H2S estimate was wrong in EVERY number. **Never
  claim a threshold crossing; never restore the eye-read nominal.**
  → `.claude/memory/prohibitions/quinone-branching.md`
- **Ascorbate on the quinone node (D-202)** — **BUILT**, wine/cascade, **default 0 and that is
  SOURCED** ("new wine has negligible ascorbic acid"); enters only via `add_ascorbate`, **never
  seed it**. **The published 2:1→1:1 signal is INEXPRESSIBLE at ANY rate** — 100× the printed
  constant still gives 0.898, because the un-dosed baseline already sits at 1.108. **Never fit
  `k_rel`.** No unsourced coefficient to relocate. Both un-built limbs push DOWN. Corpus searched:
  24 texts, 203 mentions, **no O2-limb rate**. → `.claude/memory/prohibitions/ascorbate-quinone-route.md`
- **Sotolon's pH + ethanol terms (D-205)** — **REFUSED, measured.** D-107's "the model has both
  quantities, so a real omission rather than an inexpressible one" rode 10 records as *buildable*:
  **expressible ≠ identifiable.** Two rival pH forms agree within **17 %**, so the mechanism is
  nearly free — the loose parameter is the **REFERENCE pH** (`k_sotolon_aldol` was fitted with NO pH
  term), and 3.4 vs 3.0 reports **two oxidised wines or none**. Ethanol is **BLOCKED on a Keq**
  (D-190 shape). **Pham 1995 UNREACHABLE, 2 hosts.** The emergent pH limb is a **DIFFERENT
  REACTION** and no substitute: 1.893 % where SO₂ is in excess, **0.003 %** where it decides.
  → `.claude/memory/prohibitions/sotolon-ph-and-ethanol.md`
- **Strecker methional split (D-101 → D-206)** — deriving `f_methional` from the amino-acid
  abundances is **REFUSED, measured**, and so is re-banding it: the abundances **cannot reach**
  that output (gate cancels the fraction; a 4× recompiled change moves it **0.06 %**), so it is
  the ONLY channel. D-101's "≈ 0.136" was a **mismatched pair** (true value 0.15152, 1.01 % off
  the shipped 0.15) and its "queued as its own beat" clause is **deleted**. Also: a **reused**
  `CompiledScenario` inherits the previous run's enables (**+10.3 %**) — a CONTRACT, "fixing" it
  fails 26 tests. → `.claude/memory/prohibitions/strecker-methional-split.md`
- **Beer acid-base — seven beats, the pH CURVE, the FRAME, the NITROGEN CHARGE (D-178 → D-183, D-207 →
  D-209)** — beer's pH is a **PREDICTION** and it now **AGREES**; malt phosphate REFUSED as the buffer;
  registries NEVER merged; the `CO2` slot is NEVER dissolved; acetic's producer is growth-linked, spike
  NOT modelled. **Never call the pH test a TRAJECTORY test** (D-207). A published pH is **DECARBONATED**
  — score only `degassed_ph_of_state` (D-208). **D-209 put the `N` pool on the CATION side and CLOSED the
  ~0.4 pH: `z̄` is DERIVED per ELEMENTAL N, never fit; symport is proton-NEUTRAL; t=0 is a RE-ALLOCATION;
  it rides D-179's gate. Its "day 1 is UPTAKE TIMING" is CLOSED at D-211 — see the row below.**
  → `.claude/memory/prohibitions/beer-acid-base.md`
- **Beer's growth RATE + uptake TIMING (D-211)** — `mu_max` **0.034**, band 0.031-0.040, MEASURED on
  Tyrell Fig. 4's **cell-count panel** (the third panel D-180/D-207 left uncropped); 0.098 was **2.88×**
  too fast — a **Droop→Monod** transfer, D-15's twin. **Never fit the pH course** (0.040 scores 8/8 and
  was REFUSED); EXTENT is a separate 6.3 % deviation; the `E_a` arm is measured and REFUSED. Band is
  **DRAWN** and narrowed 7.00→1.29× from the fast end. **D-210's two terms are NO LONGER PARKED —
  BOTH CLOSED at D-214** (antiport unsourced, trub pre-pitch); §9's numbers → 0.0274/0.0086. D-183 FLAGGED.
  → `.claude/memory/prohibitions/beer-growth-rate-and-uptake-timing.md`
- **Beer's WORT OXYGEN (D-213)** — **BUILT**: `o2` is seeded at **6.75 mg/L** (band 5.5-8.0, both
  edges PRINTED) and stripped by the yeast in the lag phase. **DELIBERATELY INERT — the owner
  chose it knowing** (all 3 O₂ consumers are aging-gated); it stops beer claiming t=0 anaerobiosis
  and defuses D-212 §7. Driver is biomass **PRESENT, never FORMED** (growth-coupled built + REJECTED
  on timing). **NOT a growth-Arrhenius target.** Isolability is exact at the DERIVATIVE, **not
  byte-for-byte** — proved MESH by CONVERGENCE. O₂→growth DECLINED; O₂→acetate stays REFUSED.
  → `.claude/memory/prohibitions/beer-wort-oxygen.md`
- **Trub settling + the peptide PAIR (D-214)** — **REFUSED, and D-209 §8's BOTH parked terms are
  now CLOSED.** Antiport is **ZERO hits in every beer text** (its 8 are LAB, in wine). Trub is
  **PRE-PITCH** (boil + chill, *"removed before the wort is fermented"*) and **already inside the
  1.18 control-wort anchor** — never call it an omission. A Process draining `peptide_buffer` is a
  **CHARGE VIOLATION** (pH 7.08 / 11.66), not a small acidifier; pre-anchor it is refused on
  **SHAPE** (day 7 is **3.40×** day 1, the inverse of the brief) and the window is **EMPTY by 9×**
  at the high `z̄` edge. The capacity/pKa pair is **incoherent but NOT fixed — its own beat**.
  D-211 §9's numbers → **0.0274 / 0.0086**. → `.claude/memory/prohibitions/trub-settling-and-the-peptide-pair.md`
- **Beer pH: the LITERATURE has none, and the acid COURSES were never scored (D-215)** — **NO model
  change.** Corpus null on a missing acidifier: Ca-phosphate is **BOIL** (pre-pitch, like trub),
  H⁺-ATPase is **LAB/wine only**, no beer BC-during-ferment, **Tyrell is still the ONLY pH curve**.
  Figs 9/10/14 interiors read at last: the three `Y·ΔS` timing errors **OPPOSE and nearly cancel**
  (succinic 25 pts LATE, malic 25 pts EARLY) — **no single fix helps all three**. The engine ferments
  Tyrell's wort **~2.8× too slow at day 2 (21.2 vs 59.4 %), never scored** — but flux-matching
  **OVERSHOOTS** (day 1 → 0.063 too *acidic*), so it is **NOT** established as the day-1 cause. pH
  probe **INCONCLUSIVE**: noise floor 3.2× the headroom. → `.claude/memory/prohibitions/beer-acid-course-timing.md`
- **Beer's ferment SPEED + the pitch its pH rides on (D-216)** — **REFUSED, tests only.** NOT
  D-211's doing (lag pre-dates it: 2.05×→2.81×) and NOT growth extent (**0.8 %**). `q_sugar_max`
  1.397 closes it and is **IN BAND** — §2.2's benchmark forbids it, breaking at **q≈0.6**; no
  in-band pair works and **not even removing `K_repression`** (79 % of it) reaches 0.594.
  Corner leaves **1.79×**, a LOWER bound (`E_a_uptake`'s low edge is D-19's debunked figure).
  **NEVER "correct" the 1.0 g/L pitch** — honest 0.5 makes the lag AND the pH worse; **D-211
  FLAGGED**, its 0.070 is conditional. Which anchor wins is **the owner's, unrecorded**.
  → `.claude/memory/prohibitions/beer-ferment-speed-anchor-conflict.md`
- **Beer's EARLY acetic rise (D-212)** — **REFUSED, pre-registered, BUILT NOTHING.** Never aim at
  Tyrell's day-1 **145.0**: the pH admits a **WINDOW** — **94.24-132.38** admitted by ALL 3 `z̄`
  arms (86-141 is the UNION, never quote it) — and 145 is **OUTSIDE at every arm**. The real
  shortfall is **+6 to +14**, a fifth of it, and hitting the
  acetic target would move the pH OUT of band. Refused on **MECHANISM, not identifiability** (any
  front-load factor 1.67-2.78 works at all arms). **Four candidates, FOUR different deaths.**
  An ADDITIVE source is forbidden (breaks day 7 → 182.6). D-211 §13's sweep claim CORRECTED.
  → `.claude/memory/prohibitions/beer-early-acetic-rise.md`
- **The nitrogen DOSE's charge (D-210)** — `add_dap` doses a **SALT** and BOTH ions are BUILT:
  `phosphate` (diprotic, both registries, NOT D-178's malt phosphate) and `nitrogen_charge_excess`
  (stores the EXCESS so 0.0 needs no sentinel; one slot, not a second N pool). **D-209 §8c sized
  the WRONG half** — the ammonium moves a dry wine's endpoint by **0.0** and owns the excursion;
  the phosphate owns all **−0.162 pH**. The 8 `N` inflows keep the average (88 % is an UN-DRAW).
  Native phosphate is out on the ANCHOR, not smallness. `set_ph` MASKED the guard's message.
  DAS and the arginine pool's own +1 stay unbuilt.
  → `.claude/memory/prohibitions/nitrogen-dose-salt-charge.md`

- **Beer's uptake TEMPERATURE sensitivity + Tyrell's FRAME (D-217)** — **MEASURED NULL, nothing built.** The
  corpus **cannot** re-source `E_a_uptake`: 75 hits/26 files, then **28 more** for the two-temperature-pair shape
  the first census *could not see* — all read, zero usable. **Never re-run without a NEW text**; never adopt the
  30,000 edge; **never read de Andres-Toro's −97 kJ/mol as a fit** (wrong sign; degenerate with `q_sugar_max`).
  Benchmark inertness is **EXACT (0.0000 d)**, so the lever's size **IS** the frame's distance from `T_ref`:
  0.0449 at 15 °C, **exactly 0 at 20 °C**. Tyrell's tube temp is **NOT PRINTED** — a FILL temp D-207 promoted;
  §3.2 **vindicates** it (its ⚠ is GROUNDS, not value) so **never "fix" it**; D-216's refusal STANDS.
  → `.claude/memory/prohibitions/beer-uptake-temperature-sensitivity.md`

- **Which yardstick beer's SPEED is calibrated against (D-218 → ADJUDICATED AND BUILT at D-223).** Foster's measured
  15 °C course + §2.2's window WIN and are ONE anchor; Tyrell's day-2 extract point LOSES. `q_sugar_max` 0.5 → **0.72**,
  band 0.3-1.5 → **0.634-0.818**, and it is DRAWN. One rate fitted at 15 °C lands **0.945 / 1.030 / 0.973×** at
  12 / 15 / 22 °C — 12 and 22 OUT OF SAMPLE — so D-220 §4's "level error, not a temperature-response one" was right and
  D-217 stays refused. **The 30 °C column inverted to 0.659× (model too FAST) and is NOT repaired.** The faster engine
  also exposed an ester formation half drawing acid from a pool that is 0 by construction; it is now funded.
  → `.claude/memory/prohibitions/beer-ferment-speed-anchor-conflict.md` (D-223 block at the top)
  **Everything below this line is the pre-D-223 state and is kept for the reasoning, not the numbers.**
  D-216 §11's open question is **answered on the literature and now open on a CONVERSION FACTOR.** Foster 2022
  (*Front. Microbiol.* 13:747546) is the first source carrying the whole tuple — 12.5 °P (within 0.003 SG of
  §2.2's wort), a **counted** pitch, two temps for the **same three Beer 1 strains**, same 1.010 target (Parker
  2008 sources it, so the brief's **finish line** is corroborated; its **duration** is not). Every third-party
  endpoint is at or faster than 5-7 d's fast edge. **THIS REPO SHIPS TWO PER-CELL DRY MASSES 5.6× APART** —
  18 pg/cell in the wine validation files, ~100 pg implied by `beer_generic.yaml`; D-216 §7's "~2× textbook"
  understates it (hence the **Flags**). Bracketed over both readings **and** both ends of Foster's sampling
  interval: the window survives in **1 of 8 cells**, needing the anomalous reading **and** a 72 h SAMPLE read
  as exact. At that reading the **shipped** model is **11 % off** a published endpoint, untouched.
  **Do NOT retire §2.2's window on this** — a rate satisfying Foster finishes a day early against both measured
  tails while still missing both day 2s (D-216 §4 unchanged). **Never cite Foster's pair as a temperature test**:
  the bound is cleared by the ENTIRE printed `E_a_uptake` band and first fires 1.43× out — so **D-217's refusal
  needs no re-opening**, and the near-equality of the ratio to the bare Arrhenius factor is a **CROSSING** at
  ≈58 kJ/mol, not a law. **The per-cell dry mass is SETTLED — next block.**
  → `.claude/memory/prohibitions/beer-speed-yardstick.md`

- **D-219 — a cell's dry mass is a UNIT DEFINITION, not a literature question; §2.2's window FAILS.**
  **4e-11 g/cell**: Coleman 2007 (wine's `Y_X/N`/`k'_d`/`mu_max` source) says *"assuming that each cell
  weighs 4 × 10⁻¹¹ g"* and COUNTED cells — every gram of this engine's biomass is a count × that. **Assert
  it EXACTLY; do not re-open as "which estimate".** Both old readings RETIRED: 18 pg unsourced, **~100 pg
  BACK-COMPUTED** from `pitch_gpl = 1.0` (a residual, not a mass). **Never band it from dry-yeast DOSING
  conventions.** **§2.2's 5-7 d does NOT survive** (3.63 d) and the model is **1.61× slow, not 11 % off**
  — but `q_sugar_max` NOT moved and `TYRELL_SCENARIO`'s **2.51×** excess NOT corrected (`mu_max` refit;
  6.58× day-2; pH 7/8→6/8). **2.51× and D-215's 2.8× COMPOUND, not corroborate.**
  [[feedback-a-constant-can-be-a-unit-definition]] → `.claude/memory/prohibitions/yeast-cell-dry-mass.md`

- **The SECOND measured beer course (D-220)** — **MEASURED, 3 corrections, nothing built.** Foster's
  **Supp. Fig. S1** (VECTOR ⇒ transcribed, never eye-read), in a paper mined TWICE: **D-216's "no second
  beer course" is FALSE.** §2.2's window is **MIS-TEMPERATURED, not refuted** — all 3 ale strains INSIDE
  5-7 d at **15 °C** but 2.91-3.76 d at 22 °C, so the duration is real and the **PAIR with 20 °C** is wrong
  (a BRACKET; the 20 °C interpolation is asserted NOWHERE). Engine **1.41/1.54/1.45× slow at 12/15/22 °C**
  = a **LEVEL** error ⇒ **D-217 VINDICATED**; §2.2 HIDES it. **30 °C is a CROSSING** (measured E_a collapses
  49.5→17.9, model holds 55.5→53.7), never a validation nor a reading of `E_a_uptake`. D-218's day-1-2 peak
  is **temperature-conditional** — model peaks later still. → `.claude/memory/prohibitions/beer-second-measured-course.md`

- **§2.2's BEER criterion is now 15 °C and an xfail (D-221)** — **SETTLED, no value moved.** The 5-7 d duration is REAL and
  the 20 °C it was paired with is not, so only `conditions` moved; window untouched, Foster's band sits INSIDE it. Engine
  **9.00 d**. **D-216 §6's decoupling lever is DEAD** — `E_a_uptake` was free only because the anchors sat at DIFFERENT
  temperatures; 15 °C is Tyrell's exactly, so 0.0000 d → **1.75 d against a 2.0 d window** (growth 0.2917). No magnitude
  argument survives. The conflict **INVERTS**: at 20 °C the criterion forbade a faster engine, at 15 °C it demands one
  (q 0.425-0.621 → 0.667-1.017) — a MAGNITUDE question now, never a direction one. **D-216's refusal of 1.397 SURVIVES**
  (outside across the whole `E_a` band). D-218 §3 / D-219 §5c **INVERT**: the surviving corner moves to the SETTLED 40 pg
  (5.42 d), CONDITIONAL on `E_a` ≥ 40,165 (69 % of band). → `.claude/memory/prohibitions/beer-criterion-temperature.md`

- **Beer's AROMA CALIBRATION — levels, coupling and the SINK (D-224→D-227 — SETTLED)** —
  "which ester calibration is wrong" is **ANSWERED: NEITHER**; quote the **3 %**, never 0.4 %.
  Pools ride growth EXTENT (**beer only**; wine's `E_a_esters` is derived FROM the flux form and
  wine is the control). **D-227 CLOSED D-226 §11's sink item**: `k_ester_volatil` was a **NINTH**
  drift-prone constant with `q_sugar_max` folded inside it — POOL registries can't see a SHARED
  constant. Drift 13.91 %→0.73 %, residue now `mu_max`'s TIMING; **never claim exact ester
  invariance**. **Never re-propose** the coupling, Luedeking-Piret, Wang's 23.7, or re-anchoring
  WINE's `k_isoamyl_acetate`. → `prohibitions/beer-aroma-calibration-levels.md` (D-227 block FIRST)
- **The growth CUTOFF behind beer's aroma taper (D-228)** — it was the calibration frame's
  **INOCULUM**, not the growth law: that frame pitched a flat **1.0 g/L** (2.5× a counted ale — the
  residual D-219 retired and D-222 fixed only in Tyrell's scenario). Now COUNTED (1.2e7/mL); levels
  invariant **6e-7** so no re-anchor, and the ONLY guard that can see it is the growth WINDOW
  (0.785→1.199 d). **D-226 §8 / D-227 §10's "still the thing that would fix the taper" is CORRECTED.**
  On Tyrell's counts the model is INSIDE the spread days 1-3 and peaks day 3 **when sampled DAILY**.
  Taper stays REFUSED (no ester course); the **1.55× extent overshoot is CLOSED as a nitrogen
  question at D-230** — see the next block. → `prohibitions/beer-growth-cutoff.md`
- **Beer's growth EXTENT overshoot (D-230, D-232)** — **NOT a nitrogen error; both D-222 candidates CLOSED.**
  Peyer T16 sums to **189-194 mg N/L** (assumed 200 = 1.03-1.06×) and has **no proline** — no FAN→YAN fix owed.
  Partition REFUSED: Tyrell's crop at 40 pg needs **20-26 % cell N** (band tops 0.14). Neither `yan_mgl` nor the
  **gram** is isolable from the rate fit (arms predicted 4 REDs, got 8 and 10). Residue is now **THREE-way** —
  71-92 pg cells, OR 44-56 % settled at peak, OR a different organism/medium. **D-232 REFUSED BOTH ways of
  deciding it**: the day1-vs-day3 settling profile **flips SIGN** across `mu_max`'s band (not inert — 0.348 at
  day 1, 0.003 at day 3), and the pH clock only re-measures D-209 §8's **+0.162** day-1 miss (= 0.148, 70 % of
  the envelope). **Never** harmonise onto Coleman's Y_X/N (→ ~9×). → `prohibitions/beer-growth-extent.md`
- **The per-member pH ANCHOR + the peptide pair (D-233)** — the anchor re-solve is **BUILT**, BOTH
  media: members began at 5.5062-5.7778 vs an anchor of 5.65 (beer) and 3.4208-3.5780 vs 3.50 (wine),
  now **2.3e-11**. **D-24's "y0 is held fixed" CORRECTED, not repealed** — scenario INPUTS still never
  sampled. **NEVER quote 1.287× as the band** (one parameter's own contribution; the band moves
  **1.008×**) — a PER-MEMBER error, argued from **t=0 only**. **Capacity half is BUILT at D-238**
  (rule 3, running FIRST — it feeds the anchor): every member's wort reproduces Peyer's 1.18, cost
  **0.0095 pH** day 14, CONVERGED, and **exactly 0 at t=0** (the cation solve absorbs it) — the
  mirror image of the anchor half. D-233's "moving the root-find into src makes the guard circular"
  **CORRECTED: it forbids DERIVING THE SHIPPED CONSTANT, which this does not.** Its defect pin
  stayed **GREEN** and was re-scoped, **not deleted** — do not obey that instruction, re-read it.
  Renamed **`y0_for_member`** at D-236 (rules: copper seed, then D-238's capacity).
  The census it left open is RUN at **D-234** and fully repaired by D-235/236 — see the next block.
  → `.claude/memory/prohibitions/ensemble-anchor-reanchor.md`
- **The compile-read ∩ SAMPLED census (D-234/235/236)** — **RUN, CLASSIFIED and now with NO LIVE row.
  Never re-propose it as unenumerated or either repair as open.** Predicate is compile READ ∩ the
  sampler's own set — **disjoint by construction** from D-153/156/157/159's drawability surface.
  **NEVER size it with a grep**: 26 hits, 21 invisible; I predicted 12-20 and got 32. **D-235** widened
  `StateMutation` to `(schema, state, params)`: `set_ph` 0.07896 → **1.98e-11**, `add_dap` **0.0108** pH
  (exactly 2 of 13 verbs read it), and the priced radius missed **28** test CALL sites. **D-236** made
  `reanchor_for_member` → **`y0_for_member`**, per-slot RULES: copper **16.65 % → bit-identical**, and a
  wine NAMING `copper_gpl` is left alone (**14.34 %**, honest — re-seeding would breach D-24). Fractions
  CLOSED at D-206; override knobs BY DESIGN (the MODE, D-164). **D-237: the census is WIRING-INVARIANT**
  (cascade/burst add nothing) and the copper repair is bit-identical under all three — but **D-234's
  "the cascade will WIDEN it" is WRONG in DIRECTION**: 15.38 % vs direct's 16.65 %, the mildest arm.
  **NEW LIVE row, and it is the census's COMPLEMENT, not a member**: `burst_antioxidant_initial` is
  banded **50×**, seeds a live pool, and NO Process declares it — the sampler can never draw it.
  PINNED, not repaired; the census it implies (banded + compile-read + undeclared) is NOT run.
  → `.claude/memory/prohibitions/compile-sampled-census.md`
- **The wort's free amino-acid buffering (D-239)** — **BUILT: D-209 §8's buffer-removal half is
  no longer open.** Three Asp/Glu/His SIDE CHAINS at Peyer's printed pKas, read off `N` (no slot,
  no Process, nothing on the nitrogen ledger). Charge is a **RE-PARTITION**: `z̄` 0.1772 →
  **0.23418** (width unchanged), halves cancel at wort pH to **0.1772 exactly**; lump **1.54807 →
  1.43506**, wort still 1.18. Cost **0 at t=0, +0.0023 day 1** (pool 70 % present — **CORRECTS
  D-209 §8's same-sign claim**), **−0.0202 day 7** ⇒ high `z̄` edge **0.0203 BELOW** Tyrell's
  floor, deliberately. **Day-1 miss now has NO named candidate.** Wine untouched, MEASURED 0.73 %.
  → `.claude/memory/prohibitions/wort-amino-acid-buffering.md`

- **The banded-and-NEVER-DRAWN census (D-240)** — **RUN, CLASSIFIED, PINNED; never re-propose as unenumerated.**
  D-237 §6 CLOSED: **32** compile-read/drawn-nowhere, **28** banded per scenario, 5 classes. D-234's **COMPLEMENT**;
  **OVERLAPS D-159** (11 of 28) — never call those two disjoint. `copper_typical` is in BOTH registries, correctly
  (aging-gated reader) ⇒ score per SCENARIO. **NEVER quote the unpaired 1.672× widening — NOISE**; paired **0.323**,
  and the *as-shipped* arm is **EXACTLY 0** everywhere. hidden/full: methanethiol **1.88**, dms 1.09, A420 **0.705**,
  beer **0.000**. Beer's 8 `*_typical_wort` **ABSORBED at t=0 (≤1.2e-13)**; `o2_wort_aeration_beer` worth **EXACTLY 0**
  (D-213); the 2 copper binding constants **lose a `min()`** (126-2200× excess) — pin the RATIO, never the zero.
  → `.claude/memory/prohibitions/banded-undrawn-census.md`
- **Seed drawability REPAIRED (D-241)** — `CompiledScenario.seed_reads`, **derived from the `y0` rules** so a name is
  drawable **iff** a rule re-seeds it; unioned into the DEFAULT branch only, BEFORE `exclude`. **6 of D-240's 8 live
  seeds now DRAWN** (`burst_antioxidant_initial` under `direct_burst` ONLY — D-147 zeroes the slot elsewhere, whole `y0`
  bit-identical, NOT a half-repair). **Tier half deliberately UNTOUCHED and still OPEN — never cite D-241 as settling it.**
  Paired: `dms` **2.83×**, methanethiol **2.07×**, burst **6.97×**; beer **1.000** = the null control. The 2 Coleman
  coefficients are **SUBSUMED** (the override they derive IS sampled, band 2.11× wider and containing theirs) — never
  propose a `values_for_member` hook for them. D-240's undrawn pin **forbade nothing** against this route.
  → `.claude/memory/prohibitions/seed-reads-repair.md`

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
  `d15{7,8,9}-*\`, `d16{0..8}-*\`, `d18{7,8,9}-*\`, `d19{0..9}-*\`, `d20{0..9}-*\` (**`d209-nitrogen-proton-exchange\` holds `peyer_full.txt`, the 243-page thesis extracted IN FULL after D-178 read one chapter of it**), `d170-q10-generalise\`/`d171-ordering-guards\`/`d175-ellagitannin-joint\` — incl. `d142-pulls\`+`d143-so2-binding\` (Miao **T2/3/4**), `d149-copper-refit\`
  (Nguyen **T3.1**), `d151-l16-ph\` (Carrasco-Quiroz **T1+2**), `d163-band-edges\`/`d165-wide-band\`/`d166-switch-census\`/`d167-edge-provenance\` (reusable harnesses); `_txt\carrascon-red-kinetics-2018.txt` = **Carrascón 2018 reds**.

## Not started (deferred tail; D-110's narrowing still unconfirmed by owner)
~~growth-linked excretion (D-49 opt B)~~ — **NOT a candidate: α-KB REFUSED (D-189), PYRUVATE REFUSED (D-195), and
the draining arithmetic is STRUCTURAL AND GENERAL to all three pools, so α-KG is refused BY IMPLICATION** —
inferred, not separately measured, and that is the only thing left open about it. Do not offer it as unbuilt; peptide pool; variety-specific DMSp;
yeast-autolysate spectrum; re-anchor `f_methional` (**only from LITERATURE — deriving it from the
model's own abundances is REFUSED at D-206, and so is re-banding it**); masking (cosα-blocked); D-55's stale Brett prose; **acetaldehyde
in maturation + the 0-vs-2.7 floor are NOT here any more — D-188 measured both**; ester `_eq` floors; pH factor for hexanoate/EtOAc (**sourcing-blocked: no per-pH series, and R&O's per-ester constants are isoamyl's**); `k_d2`; adduct release; closure OTR(T) (**the bottling burst is BUILT — D-187**); **residual copper is BUILT — D-191; osmotic inhibition BUILT — D-192; the post-Fenton O₂ draw is BUILT — D-196, and it rode this list for 20 records after its source was already on disk**; D-143/4 ← D-145; ~~NEW at D-202: sotolon from ascorbate via 2-ketobutyrate~~ — **CLOSED at D-203, REFUSED on identifiability** (~10 % molar conversion is its target). This line said "new/open" for 4 records after its own refusal shipped; a **STALE-LIST instance**, cf. [[feedback-check-the-blocker-is-still-blocking]]. ~~NEW at D-233: the census of parameters read at COMPILE time that are ALSO in the sampled set~~ — **CLOSED: RUN at D-234 (32 names, all classified) and both LIVE rows REPAIRED at D-235/D-236.** Do not re-open it from this line [[feedback-a-parameter-can-be-pinned-and-drawn]].

## Standing rule
- **NEVER put a whole-file line cap back (D-177, corrects D-169).** Raised 4× (150→300), then **REMOVED, not
  raised a 5th**: at 300 the file **re-pinned — exactly 300 on 8 of 9 commits — with shape GREEN everywhere**,
  because a per-block cap bounds what each record **ADDS**, and nothing bounds **how many**. D-169's licensed
  ⚠-COLLAPSE RETIREMENT **cannot pay for it** — 18 blocks / **59 lines** cite only corrected records but nearly
  all are prohibitions their corrector **SHARPENED**, so recovery is **5-15 lines** (1-4 sessions). A **DERIVED**
  total (8×blocks) is **VACUOUS** (the file runs ~3.8×). Totals are now **REPORTED**: no threshold ⇒ no target.
  A **digit-density check was REJECTED on measurement** (reads block SIZE; penalises corrected-value guardrails).
- **The per-subject detail is SPLIT OUT and stays out (D-185).** Growth was never in the live frontier
  (steady **6-17** entries citing the last 5 beats, since June) but in the **SETTLED tail** (**0 → 22**
  entries citing records 40+ beats back, monotone). That tail is **INCOMPRESSIBLE** — all 22 read live —
  so **age does not measure settledness**, and a tail budget was DESIGNED and **REJECTED**: it would
  demand deleting the longest-surviving guardrails and could never be satisfied. Fix is **granularity,
  not eviction**: ledger resident, detail by path, **NO `MEMORY.md` row** ⇒ boot unchanged (~27.8 KB).
  **Never index the `prohibitions/` files there**, and never inline a subject back into this file.
- **`CLAUDE.md` is now MEASURED too** (14 lines/block, its own shape — 8 would fire on real docs). It grew
  **66→138 lines**, **+30 the day the memory cap moved**, and now carries prohibitions. **Still open channels:
  the GLOBAL `CLAUDE.md`** (outside the repo, out of scope) and **`MEMORY.md` row COUNT** (5→40 rows, +1/record).
  **The hook cannot see growth RATE** — the +30 day arrived as two ~10-line blocks and would not fire.

**Direction is the owner's call, every time** — ask before picking the next milestone/beat, offering only UNBLOCKED options (D-66, [[feedback-discuss-disagreements]]).
