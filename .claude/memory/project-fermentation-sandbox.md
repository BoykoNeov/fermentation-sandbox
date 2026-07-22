---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status + where the canonical decision/architecture records live"
metadata:
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
  modified: 2026-07-22T07:45:02.717Z
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation
engine in Python (uv, scipy/numpy/pydantic). Repo:
https://github.com/BoykoNeov/fermentation-sandbox (default branch `main`).

**This file is session-boot context, NOT a changelog.** The per-decision
narrative lives in the repo and is the single source of truth — do NOT copy it
here (see [[feedback-batch-end-ritual]]). This file is capped at 200 lines
(raised from 150 on 2026-07-21) because it regrew 114× after the 2026-07-02
re-cut; keep entries to a name plus a one-line hook and let the reader open the
D-record.

- `docs/DECISIONS.md` — every decision D-1 … D-132, with fork, choice, reasoning.
  **The canonical archive.** Its `<!-- END INDEX -->` TOC lags — it stops at
  D-129; the entries themselves are current (body reaches D-132).
- `docs/ARCHITECTURE.md` — layering, package map, core/runtime/scenario seams.
- `docs/plans/milestone-*.md` — task checklists per milestone.
- `CLAUDE.md` — prime directives (tiers, provenance, one-directional deps).

**Status (2026-07-22 — at D-132):** M0/M1/M2 **complete**; **Milestone 3**
(sensory/OAV + Tier-3 aging, owner's pick at D-66) in progress. **D-132** made
`PhenolicBrowning`'s O2-consumption rate phenolic-driven (Ferreira 2015):
`k_browning_eff = k_browning_base + k_browning_phenolic·(tannin+anthocyanin)`,
additive not proportional, isolable at zero phenolics, first-order-in-o2
re-confirmed not replaced. **D-131** built MLF
medium-chain-fatty-acid inhibition (`mcfa` inert must-input slot, bacteriostatic
`g_FA` gate on conversion/growth NOT death, Lonvaud-Funel 1988 Table 6, plausible).
**D-130** built botrytis SO2 binding (`oxofructose` must-input, 4th carbonyl, Barbe
2000). ruff + mypy (109) + full pytest all green (1256). Collection counts drift
±1–2 between beats; not a real delta.

**DONE — do NOT re-propose these (I did, twice, from stale "Next:" breadcrumbs;
see [[feedback-verify-latest-state-not-breadcrumbs]]):** ALL lumps are speciated —
esters→3 (D-96), **fusels→5** (D-99: propanol/isobutanol/active_amyl_alcohol/
isoamyl_alcohol/2_phenylethanol), amino_acids→8 (D-100), mercaptans=methanethiol
false-lump (D-110). **"No pool in the project is lumped any more" (D-110).**
**Beat 1b (descriptor projection) is COMPLETE** — slice 1 (D-95, `descriptors.py`,
max rule) + slice 2 (D-98, `compression.py`, Stevens + flip-sensitivity); D-98
*retired* the "Next: beat 1b slice 2" item, its only remaining piece (masking) is
`cosα`-blocked. Milestone-3's substantive builds are largely done; **most
remaining work is blocked on external sourcing**, not on more building.

**Live threads** — full reasoning in the named D-record, not here:

- **Phenylalanine / 2-PE de-novo route — SHIPPED (D-118), then MEASURED (D-119).**
  `f_non_ehrlich_phenylalanine` = 0.975, `f_de_novo_2_phenylethanol` = 0.9827.
  D-119 read Minebois Fig. 6A and **changed nothing**: it corroborated the 2.5%
  numerator in µM independently of the algebra, but the in-study 0.963 has a
  *different denominator* (her 109 µM total 2-PE vs the model's Wang-anchored
  ~235 µM). **The blocker moved, it did not lift** — from "figure unread" to
  "does the Phe flux scale with total 2-PE?". Residual risk is **guard-safe**.
  **Never put 0.963 in a sampled field.**
- **Three D-119 errors worth not repeating.** (1) A de-novo fraction is **not
  scale-invariant** — denominators don't transplant between fermentations.
  (2) The "it's only a T4 snapshot that climbs later" defence is **refuted**;
  the in-study fraction is flat. (3) **Fig. 6A's bar semantics are MIXED** — for
  2-PE the large number is the TOTAL, for isoamyl/isobutanol the UNLABELLED
  segment — and the first pass asserted a *uniform* rule, then wrote
  "2.57% ≈ 2.5%" to protect it. **Sanding a rounding mismatch to keep a tidy
  rule is the tell**; derive each bar from its own printed constraints. The
  de-novo share is `1 − I.E` and robust to the split either way.
- **Ask Querol** (`aquerol@iata.csic.es`) for raw SI: Phe dose vs total 2-PE
  across conditions. DataSets are "upon reasonable request" (verified by
  fetching). Isoleucine's ¹³C tracer rides along on the same request.
- **Ramey & Ough 1980 was READABLE after all (D-123) — the ACS paywall was not
  the paper's only door.** It is an open scanned PDF on the author's winery site
  (`rameywine.com/wp-content/uploads/2023/04/ester_hydrolysis.pdf`); read by
  rendering PDF pages. D-121 recorded it "blocked / not a search problem" against
  the ACS `CLOSED` copy — a fresh *secondary/author-hosted* search found it.
  **Lesson (D-118/D-119 inverted): don't record one host's paywall as the whole
  paper's blocker — check author-hosted/secondary copies first.** Still
  paywalled: — (Makhotkina & Kilmartin 2012, PMID 22868118, is now FETCHED and
  spent at D-126; see the ester-aging bullet).
- **Single-host obligation still OPEN.** Everything Minebois rests on the PMC
  deposit; Wiley full text and PDF are HTTP 402. The figure image is the *same*
  deposit, not a second host. Two live parameters ride on that one figure.
- **Isoamyl's de-novo entry — REFUSED at D-120, measured not built. Do not
  re-open it as a build.** The cap is `gate *= (1 − f_de_novo)`, a
  **one-directional ceiling**, so it only helps where the model OVER-attributes.
  Every alcohol is on the wrong side (model vs Minebois AA share: 2-PE 0.87/3.72,
  isoamyl 2.94/5.34, isobutanol 5.14/8.78). **D-118's "class of error" claim is
  overturned in the actionable direction** — all three *are* de-novo dominated
  in-study, but the model is *more* de-novo than measured, not less. Two further
  kills: the cap can't reach isoamyl's valine branch (primary loop only, and it's
  the larger branch); and it's a **rate knob on a supply-limited quantity** —
  every precursor exhausts, so the share sits on D-112's `(1−f)` ceiling.
- **The shipped 2-PE cap is INERT at the realistic dose** (0.871% either way) —
  it bites only at 4× dose where Phe stops exhausting. D-118 is **not** wrong:
  the route's measured job is the instantaneous **carbon-refund guard**
  (1.125×→0.584×, pinned by a counterfactual test), a different quantity from the
  integrated share. What doesn't survive is reading the route as the fix for the
  18.9% over-attribution — at this dose `f_non_ehrlich_phenylalanine` 0.53→0.975
  did that. **Never conflate the two justifications again.**
- **The live lever on isoamyl is `f_non_ehrlich_leucine`, not a route** —
  Crépin's 0.815 (shipped) vs the ~0.29 implied by Minebois (the open D-103
  conflict). Her `f` and her de-novo share are **one study**: their consistency
  is internal, **not** corroboration.
- **Ester aging — three ester Processes now ship (full detail in DECISIONS.md D-123→D-127).**
  (1) `EsterHydrolysis` isoamyl-acetate k/Ea re-anchored to Ramey & Ough 1980 real wine (D-123),
  then pH-explicit: first-order [H⁺] (D-124) → full multi-species tartrate law (D-125,
  `h=N(pH,tart)/N(pH_ref,7.5)`, `k_H2T` NEGATIVE shipped faithfully, wine-only). (2)
  `EthylHexanoateHydrolysis` (D-126), Makhotkina & Kilmartin 2012 real-wine k_obs+Ea=68kJ/mol,
  floored+grafted, NO pH factor (deferred). (3) `EthylAcetateEsterification` (D-127), bidirectional,
  model-derived. All SPECULATIVE (D-1 floor). Deferred: composition-dependent `_eq` floors (all
  three); pH factor for hexanoate/EtOAc (R&O Table VII ratios are isoamyl's, not ported). Direction
  facts (still true): MCFA esters hydrolyse WITH the acetates; branched/polyprotic ethyl esters
  (diethyl succinate, ethyl lactate) are the forming family — sim tracks none. Ethyl decanoate trap:
  fast decrease, no uniform C4–C10 constant.
- **D-127 book/paper sweep — findings durable at `M:\claud_projects\temp\ferment\_findings\`**
  (book-sweep.md + obtained-articles.md). Beyond the EtOAc build, it stocked **anchors for future
  axes — don't re-fetch:** oxidation O2-consumption (Ferreira 2015) STILL UNBUILT. **MLF fatty-acid
  inhibition (Lonvaud-Funel 1988) = D-131 SPENT it** — `mcfa` inert must-input + bacteriostatic
  `g_FA` gate; book-sweep note was WRONG (C10 23µM = "marked delay" not total; C12 white = total),
  primary Table 6 is the calibration; B (yeast MCFA synthesis) deferred (per-yeast yields = D-121).
  Findings: `_findings/D-131-mlf-fatty-acid-xval.md`. **SO2 carbonyl binding = D-128 pattern AGAIN (~80% built at D-51);
  D-130 SPENT it** — Barbe's real add = botrytis carbonyls as a must INPUT (`oxofructose` inert slot,
  4th carbonyl, plausible); DHA omitted (transient), gluconolactone deferred (D-121).
  Beer 3-sugar (glucose/maltose/maltotriose) kinetics NOT present in the 5 beer books. Rayne 2016's
  CALC ester k runs 6–18× off R&O MEASURED — never an absolute.
- **Ethanol ceiling — D-128 SPENT Luong (provenance re-anchor only); D-129 BUILT the gap. Don't
  re-propose either.** Coleman death (D-13) is linear/unbounded (valid ≤300 g/L) → high-sugar musts
  over-fermented. **D-129 `EthanolToleranceDeath`** adds `Φ=k_d2·max(E−E_tol,0)²` INSIDE the death
  mechanism (D-13-compatible, not the retired uptake wall) → musts STICK sweet; EXACTLY zero below tol
  → in-envelope byte-for-byte. WHERE=sourced 142, HOW-SHARP=spec `k_d2`. Deferred (in DECISIONS):
  osmotic inhibition >~200 g/L; measured `k_d2`.
- **D-104's inverted split / un-inversion build** — scoped, UNSOURCED, not
  started, owner's call. D-116 moved its gate off the transaminase rate onto
  **in-situ [E] + de-novo-KIC and decarboxylase fluxes**; open by trying to
  source those three, and by pricing the **Minebois/Crépin leucine conflict**
  (29.3% vs 77–86% to protein — two bands, never averaged, D-103).
- **Closed, do not reopen:** the leucine shortfall (D-112, measured not built);
  shared-BAT parsimony (D-116, a parsimony *loss*); Rollero's isoamyl-acetate
  enrichment (D-115, built — ratio 1.05 is the deliverable); no pool is lumped
  (D-110, `lumped` kept **dormant**, not deleted).
- **Recorded, not tuned:** realised Phe share under-shoots (guard-safe); static
  share ignores feedback inhibition (barely reachable); de-novo decarboxylation
  CO₂ uncharged (widens a D-19 gap, ~3e-5, both media); ester/alcohol ratio
  marginally above 1.
- **Deferred tail** (D-110 narrowed it and shipped one item; the narrowing is
  **still unconfirmed by the owner**, rest NOT started): `oav`→`magnitude`;
  Pham's pH + ethanol terms; sotolon enantiomers; methionine sink + `methionol`;
  growth-linked excretion shape (D-49 option B); peptide pool; variety-specific
  DMSp; sourced yeast-autolysate spectrum; re-anchor `f_methional`; masking
  (blocked on cosα); D-55's stale Brett-phenol prose in `chemistry.py`. Plus
  closure O₂ ingress (load-bearing on the sotolon headline) and acetaldehyde in
  maturation + the 0-vs-2.7 unsulfited floor.

- **D-132 SHIPPED — Ferreira 2015 phenolic-driven O₂ rate. Do NOT re-propose.**
  `k_browning` renamed `k_browning_base` (unchanged, 3.0e-4/h) + new
  `k_browning_phenolic` (1.1e-3 L/(g·h)) boost on `tannin+anthocyanin` (D-79,
  guarded absent on beer). First-order-in-`[o2]` re-examined at the 2026-07-22
  revisit and RE-CONFIRMED (Ferreira's R²>0.989 "linear" headline is a
  cross-cycle re-saturation artifact; the paper's own within-cycle Ln[O₂]-vs-time
  diagnostic is first-order) — no MM/`Km`. Additive baseline+boost (never pure
  proportionality — else white/beer browning → 0). Boost is browning-side only,
  NOT `k_ethanol_oxidation` (documented caveat: acetaldehyde may be
  under-produced at high phenolic load). Typical red (2.3 g/L tannin+
  anthocyanin) now lands ~0.58 mg/L/day at fresh 8 mg/L O₂ — in Ferreira's
  0.5–0.7 band (was ~0.1, a ~6–8× error). Isolable at zero phenolics (GATE-1).
  Full receipts + 5 new tests in `docs/DECISIONS.md` D-132.
- **NEXT UP — D-133 initial-burst antioxidant pool, NOT YET BUILT.** Finite
  non-SO2 antioxidant pool for Ferreira's day-1 fast rate (0.54–8.2 mg/L/day,
  R²=0 vs average, Cu-positive/phenolic-negative) — OWN pool, must NOT reuse
  D-132's phenolic driver and must NOT double-count `SulfiteOxidation` (D-72).
  D-128 re-anchor pattern; all O₂ outputs stay SPECULATIVE. Design handoff still
  at `M:\claud_projects\temp\ferment\_findings\D-132-133-ferreira-o2-consumption-design.md`.
- **Standing rule that outlived its bullet (D-66):** direction is the owner's call,
every time — ask before picking the next milestone (see
[[feedback-discuss-disagreements]]).
