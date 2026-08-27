---
name: wine-nitrogen-budget
description: "D-243/D-244 — wine's nitrogen budget is AUDITED and its two channels now PARTITION: yan_mgl is the total, the yield fit is evaluated there and held at Coleman's edge, and six fusel guards are strict xfails naming the cost"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — wine's nitrogen budget (D-243, REPAIRED at D-244).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about wine YAN, `yan_mgl`, `amino_acids_gpl`, `biomass_N_fraction`,
Coleman's `Y_X/N` regression, the fusel node's de-novo shares, or "the nitrogen budget is wrong".
Every bullet is *what it forbids* + the record to read for *why*. **If a prohibition looks
unconvincing, go read D-243/D-244 — do not argue past it from here.** Siblings:
`seed-reads-repair.md` (D-241), `banded-undrawn-census.md` (D-240).

- **The evaluation point is CLOSED — never re-propose it as open.** D-243 found it and left it on
  the owner's call; D-244 took the call. `yan_mgl` is the must's **TOTAL assimilable nitrogen** and
  the eight amino-acid pools are **carved out of it**, implementing D-32's own premise. Both
  records' "OPEN, owner's call" text is spent. `Corrects:` D-14 and D-32.

- **NEVER re-propose "sum the channels into the fit".** It does not merely extrapolate — it does
  not RUN. `biomass_N_fraction` is built as a banded `Parameter` and is **refused above a total of
  444.0 mg N/L**, so 35 suite scenarios stop compiling. The only way through is widening the
  `[0.03, 0.15]` bracket, which `seed-reads-repair.md` forbids: it is what the ensemble draws for
  the constant governing biomass. It was also rejected on principle — D-243's finding is a
  *declaration* defect, and summing leaves `yan_mgl` named for a total it does not hold.

- **The high hold is EPISTEMIC. Never restate it as saturation.** Above
  `biomass_N_yield_fit_yan_max` (350 mg N/L, Coleman's own top treatment) the fit is evaluated at
  the edge and the derived tier drops to SPECULATIVE. No source here says `Y_X/N` plateaus; a first
  draft justified the hold by an invented physiological ceiling on cell nitrogen and that is the
  D-203/205/206 error. A sourced high-nitrogen yield curve would REPLACE the hold, not confirm it.

- **The LOW edge is recorded and deliberately NOT enforced — do not "finish" it.** `f_N` is
  monotone in YAN, so below the span it only falls, infimum `1/exp(3.50) = 0.0302`, inside the
  bracket and inside physiology. A low hold would also move **Varela's 50 mg N/L arm** — the
  project's only independent wine dataset, D-56's firewall. Guarded in both directions.

- **There is NO minimum-ammonium floor and adding one needs a source.** Real musts always carry
  some ammonium; nothing states how much. An ammonium-poor must stays legal and its residual sugar
  is the honest output. The refusal line is exactly `amino-acid N > yan_mgl`, never a clamp to zero.

- **Migrate a pre-D-244 scenario by ADDING its dose's nitrogen to its declared YAN** —
  `amino_acid_dose_nitrogen_mgl` computes it and the pitch state comes out bit-identical. Do NOT
  re-author a fixture's composition instead; a draft did, and produced a headline resting on a dose
  nobody published [[feedback-migrate-a-fixture-by-its-state-not-its-intent]].

- **The six fusel xfails are STRICT and the 80 % floor is NOT to be lowered.** De-novo synthesis is
  growth-linked; the corrected yield roughly halves biomass, so propanol falls to 77.4 % and
  isobutanol to 76.0 % against Crépin/Rollero's sourced floor, and isoamyl attributes 5.42 % to
  amino acids against Minebois's 5.34 %. Closing that re-opens D-109's supply premise and the D-120
  no-cap refusal — a fusel-node beat, on the owner's call, never absorbed into another.

- **Both source-commensurate fixtures were violating their own comments, and that is RECORDED not
  repaired.** Crépin's probe ran at 405.4 mg N/L against a 180 mg N/L paper (2.25x) and Rollero's at
  475.4 against 250 (1.9x), for the whole of D-109→D-115. Fixing it needs those papers' synthetic-
  medium compositions; substituting a grape-must partition is a category error
  [[feedback-a-generic-partition-is-not-a-defined-medium]].

- **Re-recorded values are NOT free to re-tolerance.** `phenylacetaldehyde`'s 2 y pin moved **97x**
  (2.30e-09 → 2.24e-07) — fewer cells drain less phenylalanine, so the aging Strecker route's pool
  survives fermentation. `faded_anthocyanin` and `ellagitannin` did not move and are the controls.
  The D-14 identity's residue grew 0.6 % → 7.8 %, autolysis leaves 49.7 % (was 45.9 %), the dry
  sotolon arms rose 1.38-1.43x while the sweet one moved 1.005x, and **Rollero's isoamyl enrichment
  ENTERED D-111's validated band** (0.018 → 0.0264). Never widen a band to swallow both readings.

- **Still true from D-243, unchanged:** the numerator is SOUND (Varela print their own proline
  subtraction, 380/65 total → 300/50 assimilable) — do not re-check it. **NEVER write the beer
  symmetry**: wine's demanded f_N sits between the engine's own two values, so a nitrogen-budget
  explanation is ADMISSIBLE here, the OPPOSITE of D-230. It is **not** a conservation defect (the
  ledger closes to 3e-14 through all 12 N-writing sites). The slope disagreement is a
  RE-EXPRESSION of D-56 finding 3, not a finding. **Do NOT tune anything against Varela, ever.**
  D-241 §2/§3's SUBSUMED verdict stands with its scope corrected to YAN 66.0-324.8.

- **Untouched and named, not oversights:** the tier half (the pools' provenance tiers still do not
  propagate into the biomass they become — D-244 adds a third instance, closes none), and
  `nitrogen_uptake_charge_wine`'s **proline tension** — its uncertainty note derives an
  ammonium-to-amino ratio from the Handbook's raw 3-10 % / 25-30 % with no proline subtraction while
  its own `conditions` field excludes proline. Pre-existing; found here, left where it was found.

Measurements: `M:\claud_projects\temp\ferment\d243-wine-nitrogen-audit\` and
`d244-wine-n-evaluation-point\` — `FINDINGS-pre.md`, `probe_band.py`, `probe_ceiling.py`,
`probe_partition_stall.py`, `probe_denovo2.py`, `probe_pins.py`, and the `patch*.py` migrations.
