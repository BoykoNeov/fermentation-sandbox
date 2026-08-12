---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-08-12T13:14:53.163Z
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation engine in Python (uv, scipy/numpy/pydantic). Repo: https://github.com/BoykoNeov/fermentation-sandbox (branch `main`).

**Session-boot context: PROHIBITIONS and POINTERS only** — not a changelog. Every bullet is *what it forbids* +
the D-record to read for *why*. If a prohibition looks unconvincing, **go read its D-record — do not argue past
it from this file.** **Per-subject detail lives in `.claude/memory/prohibitions/` and is reached BY PATH from the
ledger below — those files carry NO `MEMORY.md` row, so they cost nothing until read (D-185).**
**Caps: 8 lines per BLOCK here, 14 per block in `CLAUDE.md`, 320 chars per `MEMORY.md`
index row — and NO whole-file total, removed at D-177** (`.claude/hooks/check_memory_size.py`;
[[feedback-batch-end-ritual]]). Distil the NEW block; **never evict an old prohibition to buy a line.**

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

## Status (2026-08-12)
M0/M1/M2 **complete**. **Milestone 3** (sensory/OAV + Tier-3 aging, owner's pick at D-66) in progress, at **D-197**;
sensory 1a/1b closed, **D-139's leftovers ALL closed** (D-148/D-149). Suite **1749**. Some
"blocked on external sourcing" items are **NOT** — D-196's source was on disk for 20 records. **Slot/Process/oxidative-set counts live in `docs/ARCHITECTURE.md` —
do NOT restate them here: that duplication is exactly what rotted the doc (D-184).**
**Beer acid-base = SIX beats (D-178 solver → D-183 acetic's rate law); D-180's BOTH omitted terms BUILT.**
**`ACID_STATE` is NO LONGER medium-agnostic** (D-179): it is *wine's* registry, beside `BEER_ACIDS`, keyed off
`StateSchema.medium`. **Beer's pH is a PREDICTION.** **Next beat is the owner's call.**

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
  state); `K`/`n` are a **derived pair, not two bands**. `_SWEET_BRIX=70` is OPEN, owner's call.
  → `.claude/memory/prohibitions/osmotic-high-sugar.md`
- **Post-Fenton second O₂ (D-196)** — **BUILT** (cascade-only, still non-default): the radical
  takes an O₂, `share × activation`. **NEVER "unsourced"** — 20 records said so, the book was on
  disk. Per-ACETALDEHYDE not per-H₂O₂; 1.0 is an **upper bound** (hydroperoxyl limb named, not
  built); Gate 1 not breached. Double-counting **died on supply-limitation**, D-141's budget
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
- **Beer acid-base — all six beats (D-178 → D-183)** — beat COMPLETE, beer's pH is a **PREDICTION**;
  malt phosphate REFUSED as the buffer; registries NEVER merged; the `CO2` slot is NEVER dissolved;
  acetic's producer is growth-linked and the spike is NOT modelled.
  → `.claude/memory/prohibitions/beer-acid-base.md`

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
  `d15{7,8,9}-*\`, `d16{0..8}-*\`, `d18{7,8,9}-*\`, `d170-q10-generalise\`/`d171-ordering-guards\`/`d175-ellagitannin-joint\` — incl. `d142-pulls\`+`d143-so2-binding\` (Miao **T2/3/4**), `d149-copper-refit\`
  (Nguyen **T3.1**), `d151-l16-ph\` (Carrasco-Quiroz **T1+2**), `d163-band-edges\`/`d165-wide-band\`/`d166-switch-census\`/`d167-edge-provenance\` (reusable harnesses); `_txt\carrascon-red-kinetics-2018.txt` = **Carrascón 2018 reds**.

## Not started (deferred tail; D-110's narrowing still unconfirmed by owner)
Pham's pH + ethanol terms; growth-linked excretion (D-49 opt B) — **α-KG ONLY and INFERRED, never
measured; α-KB REFUSED at D-189 and PYRUVATE at D-195, both spent**; peptide pool; variety-specific DMSp;
yeast-autolysate spectrum; re-anchor `f_methional`; masking (cosα-blocked); D-55's stale Brett prose; **acetaldehyde
in maturation + the 0-vs-2.7 floor are NOT here any more — D-188 measured both**; ester `_eq` floors; pH factor for hexanoate/EtOAc (**sourcing-blocked: no per-pH series, and R&O's per-ester constants are isoamyl's**); `k_d2`; adduct release; closure OTR(T) (**the bottling burst is BUILT — D-187**); **residual copper is BUILT — D-191; osmotic inhibition BUILT — D-192; the post-Fenton O₂ draw is BUILT — D-196, and it rode this list for 20 records after its source was already on disk**; D-143/4 ← D-145.

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
