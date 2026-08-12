---
name: residual-copper-fining
description: "Residual copper after fining (D-191) — the credit is BUILT, the sulfide half is measured and NOT built"
metadata: 
  node_type: memory
  type: project
  originSessionId: 42d75f5f-6caf-4f0a-904b-2781dd61f15c
  modified: 2026-08-12T07:18:03.101Z
---

**Live prohibitions — residual copper after fining.** Split out at D-185's pattern; the status
ledger points here by path. Read it when working on this subject. Every bullet is *what it
forbids* + the D-record to read for *why*. If a prohibition looks unconvincing, **go read its
D-record — do not argue past it from this file.** **Never evict an old prohibition to buy a line.**

**The copper credit is BUILT (D-191). `add_copper` writes the `copper` slot. Do not re-propose
it as unbuilt, and do not "restore" the removal-only behaviour.**
- **D-44's "the CuS drops out with the lees" is FALSE and CORRECTED** — UWC 2nd ed. Ch. 24 names
  it an incorrect older-textbook assumption; the products are dispersed Cu(I)-sulfhydryl
  nanoparticles. **D-149's "the two coppers never meet" is CLOSED.** The verb docstring's
  *"nothing sources a residual-copper fraction"* was **false when written** — the source was
  already in the project's own library.
- **`copper_fining_residual_fraction` = 0.95, band [0.95, 1.0]. NEVER midpoint it to 0.975**
  (invents a central estimate no source contains) and never widen it: 0.95 is **PRINTED**, 1.0 is
  **CONSTRUCTED** (a dose cannot be more than fully retained). The whole band is worth **1.5
  points of a ~30-point effect** — the coupling is the finding, the nominal's position is not.
  Nominal sits ON the low edge, so a band-edge screen finds that arm **bitwise identical** to
  nominal: arithmetic, not inertness. It is in `test_switch_site_census`'s exact set **with all
  three `add_copper` parameters**, each there for one shared reason — one sourced edge, one
  constructed edge, no interior any source vouches for.
- **Clark et al. 2015 is RENDERED, NOT READ** — the >95 % is UWC's rendering of it. Never cite it
  as agreeing or disagreeing beyond the sentence UWC prints (the D-190 amendment's rule).
- **The sulfide half is MEASURED and NOT BUILT — do not call it a caveat and do not build it
  without the owner.** Routing fined H₂S into `bound_h2s` nearly **doubles** the reservoir
  (19.70 → 39.19 µg/L) but moves free H₂S only **25.415 → 26.505 µg/L (+1.09, 1.043×) at 3 y**,
  because release runs 1.9 %/yr. **`BoundH2SRelease` is "not copper-coupled, on purpose" (D-135)**
  — but that refusal rested on a PLS coefficient not being a stoichiometry, and the fining path
  **has** one, so this is a LIVE candidate, not a closed one.
- **D-45's mercaptide carbon flow is FLAGGED, not fixed** — it books carbon as *leaving* the wine
  on the same retracted precipitation mechanism; honest destination is `bound_methanethiol`.
- **k_copper_multiplier's D-159 freeze was SCENARIO-bound and is now THAWED** — 0.0 unfined,
  6.62e-4 fined, 1.27e-2 fined+O₂ (control `mu_max` 42.4). D-159's prose diagnosis was **right**
  and is now obsolete; its two docstrings went **stale, not red** (no assert). **Nothing here
  re-opens `k_copper_multiplier` (§2.5 CLOSED, D-149) or its D-154 band** — this moves the slot
  the multiplier READS.
- **A fixed O₂ dose BOUNDS total browning** — fined:unfined A420 runs **1.267 early → 1.030 at
  the end**. The small end-of-run gap is the browning route winning a bigger O₂ *share*, not a
  bug: **do not "fix" it.** Both ends are pinned.
- **No legal gate, deliberately** — 0.5 mg/L lands under the EU 1 mg/L limit, larger doses over
  it; over-fining is real practice and the sim is not a regulator. Handbook of Enology prints the
  CuSO₄ allowance as **2 g/hL (§4.7.1) and 1 g/hL (§8.6)** — an inconsistency inside one book;
  cite both or neither.
- **Three named extrapolations, none measured:** fraction applied to the whole dose (matches the
  source's own design — its wines were dosed *with* equimolar H₂S); white wine → red; `f_copper`
  linear past the natural-wine copper spread. Plus an **activity caveat the model cannot
  express**: complexed Cu is "non-labile" until O₂/Fe(III) re-oxidises it, and the slot is TOTAL
  copper (the basis `k_copper_multiplier` was calibrated on).

Receipts: `M:\claud_projects\temp\ferment\d191-residual-copper\` — `PREREGISTER.md`, `RESULTS.md`,
`probe1_browning.py`, `probe2_sulfide_half.py`, `probe3_thaw.py`.
Related: [[feedback-a-scope-note-can-carry-a-mechanism-claim]],
[[feedback-check-the-blocker-is-still-blocking]],
[[feedback-a-null-result-needs-a-positive-control]],
[[feedback-nominal-on-a-band-edge-is-not-inertness]].
