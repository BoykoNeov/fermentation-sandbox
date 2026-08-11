---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-08-11T22:10:25.435Z
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
M0/M1/M2 **complete**. **Milestone 3** (sensory/OAV + Tier-3 aging, owner's pick at D-66) in progress, at
**D-185**; sensory 1a/1b closed, **D-139's leftovers ALL closed** (D-148/D-149). Suite **1609**. Most remaining
work is **blocked on external sourcing**. **Slot/Process/oxidative-set counts live in `docs/ARCHITECTURE.md` —
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
  reparameterised; oak's share is a JOINT. **ZERO edges moved** bar `E_a_oak_extraction`.
  → `.claude/memory/prohibitions/band-edges-and-provenance.md`
- **Oxidation — direct set + burst (D-132 → D-137, D-147, D-149 → D-152)** — burst is **WIRED and
  NON-DEFAULT**, not dead code and not refused; copper 600 L/g CLOSED and BOUNDED; a pH term on
  activation is REFUSED; the O2 gate is Fe(II)+O2 with SO2 INVERTED in role.
  → `.claude/memory/prohibitions/oxidation-direct-and-burst.md`
- **Oxidation — cascade + SO2 dosing (D-141 → D-148)** — cascade is **BUILT and NON-DEFAULT**, benchmark
  ACTIVE not open work; must-dosing REFUSED; "the sim under-binds SO2" WITHDRAWN; §2.4 CLOSED.
  → `.claude/memory/prohibitions/oxidation-cascade-and-dosing.md`
- **Aroma compounds + Milestone-3 tail (D-117 → D-136, D-146, D-176)** — esters/fusels/2-PE/sulfides/
  closures; methionine sink is **BLOCKED not deferred**; the `oav`→`magnitude` rename must NOT be
  "finished"; EtOAc eq is Berthelot-coupled. → `.claude/memory/prohibitions/aroma-and-milestone-3-tail.md`
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
  `d15{7,8,9}-*\`, `d16{0..8}-*\`, `d170-q10-generalise\`/`d171-ordering-guards\`/`d175-ellagitannin-joint\` — incl. `d142-pulls\`+`d143-so2-binding\` (Miao **T2/3/4**), `d149-copper-refit\`
  (Nguyen **T3.1**), `d151-l16-ph\` (Carrasco-Quiroz **T1+2**), `d163-band-edges\`/`d165-wide-band\`/`d166-switch-census\`/`d167-edge-provenance\` (reusable harnesses); `_txt\carrascon-red-kinetics-2018.txt` = **Carrascón 2018 reds**.

## Not started (deferred tail; D-110's narrowing still unconfirmed by owner)
Pham's pH + ethanol terms; growth-linked excretion (D-49 opt B); peptide pool; variety-specific DMSp;
yeast-autolysate spectrum; re-anchor `f_methional`; masking (cosα-blocked); D-55's stale Brett prose; acetaldehyde
in maturation + the 0-vs-2.7 floor; ester `_eq` floors; pH factor for hexanoate/EtOAc; osmotic inhibition >~200 g/L; `k_d2`; adduct release; closure OTR(T) + bottling burst; no post-Fenton O₂ draw (D-142); **`add_copper` never writes the `copper` slot** (needs a residual-Cu fraction); D-143/4 ← D-145.

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
