---
name: wine-nitrogen-budget
description: "D-243->D-248 — wine's nitrogen budget: the two channels PARTITION (D-244), uptake is UN-COUPLED from growth demand (D-248, four fusel xfails closed), and yan_mgl's two meanings are measured-and-REFUSED"
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

- **The fusel xfails are now FIVE (D-245 §8 corrects its own §6's "four"), and the 80 % floor is still NOT to be lowered (D-245 measured
  what D-244 asserted).** One knob reproduces every crossing — hold the fixture and sweep only
  `biomass_N_fraction` 0.0578 → 0.1068 — so the must's 2.25x nitrogen is a CONSTANT, not a cause.
  **The leg that moved is the DENOMINATOR**: threonine/valine/leucine consumption is identical to
  5 dp at every biomass (those pools exhaust either way), while the alcohols halve. The same fix
  removed a **~2x over-production** (isoamyl 353.5 → 191.5 mg/L against a 172 anchor), so the floor
  had been passing on inflated de-novo carbon. **Isobutanol's miss was a HARNESS BUG, not a model
  gap** — the D-109 helper charged it valine carbon that becomes isoamyl (2.533x, live since
  D-111); measured on the residue the model routes it is 90.5 % de novo and its guard is GREEN.
  Propanol's 77.4 % is the real miss; closing it needs Crépin's medium.

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

- **The availability gate's commensurability defect (D-246 §6 → D-247) — MEASURED and REFUSED, and
  D-247 CORRECTS the record that proposed it.** `depletion_gate` scales its half-saturation by the
  must-SPECTRUM share, and a per-species override breaks the premise that cancellation rests on —
  a real defect. But the **composition-only** correction (each share re-referenced to the declared
  must, `Σf` preserved so the level cannot move) is worth **0.796275 → 0.796017** on propanol:
  6.9 % of the gap to the 0.80 floor, **in the wrong direction**, stable across rtol 1e-6→1e-10,
  and ≤0.5 % on every fusel with all Processes live. D-246's probe cleared the floor only because
  scaling `K_amino_acids` uniformly is a **LEVEL** change — the gate reads only the product
  `K·f_i`, so it is bit-identical to scaling all eight shares, and it varied nothing *relative*.
  **The residual propanol miss is therefore UNATTRIBUTED, not gate-attributed.** The gate's
  reference is an unsourced modelling device either way — D-100 argued the spectrum scaling as
  dynamic range with zero new parameters and declined per-species Michaelis constants as the D-98
  trap — so **only literature half-saturations settle it: a sourcing ask, never a core change.**
  Do not re-propose the repair. If it is ever revisited, D-247 §5 prices the compile seam: mint
  gate-scale parameters, **never overwrite `must_aa_fraction_*`** (`spectrum_carbon_per_nitrogen`
  reads those same shares and D-34's structurally-non-negative debris carbon leans on the ratio),
  conditional on all eight pools being overridden so a "hold the must, spike the leucine" scenario
  is untouched — and note that re-referencing erases, by construction, the very effect a spike
  scenario exists to study. [[feedback-a-uniform-rescale-cannot-test-a-composition-claim]]

- **Assimilable-N uptake is UN-COUPLED from growth demand — BUILT at D-248. Never re-propose the
  40.8 % residual, and never re-propose it as a parameter.** The cause was arithmetic: the D-32
  swap ran at `psi*gate*f_N*base_dx`, strictly **below** growth's own `f_N*base_dx` draw, so
  ammonium could only fall — and at zero, growth's Monod stopped growth and the swap (proportional
  to `base_dx`) stopped with it, freezing the pools at 61.394 % of initial.
  `AssimilableNitrogenUptake` draws the identity-agnostic pair at `r*mu_max*f_N*X*gate(aa)` — a
  **capacity**, reading `mu_max` as a constant and **never** `biomass_growth_rate`. Residual
  **40.8 % -> 0.62 %** (Crepin 0.2 %).
  - **The load-bearing anchor is INTERNAL and unfitted**: the seam already sets
    `f_N = 1/Y_X/N(N_init)` so biomass lands at `Y_X/N x N_init`, an identity needing complete
    consumption. **61.6 % -> 98.4 %.** Crepin's 0.2 % is the independent check, not the target.
  - **The shipped `r = 1.0` is a BOUND, not a level.** Residual, biomass and every fusel share are
    unmoved across a **200x** sweep; the residual is set by `K_amino_acids`'s asymptote, not by
    `r`. **The TIMING is separate and it MISSES** — 18.6 h to 90 % consumed vs Crepin's N_T of
    28 h. Never read "insensitive across 200x" as "the time course was checked".
  - **The carbon must NOT go to sugar and this is not stylistic.** The swap's no-hexose guarantee
    is a property of its RATE being proportional to growth's draw; an un-coupled flux has no such
    bound and would create hexose at `base_dx = 0`. Skeleton parks in `amino_acid_skeleton_carbon`
    — **elemental carbon, weight 1.0** (the `N`-slot idiom, NOT the `debris`/glucan one).
  - **FOUR strict fusel xfails CLOSED and no sourced threshold moved.** Propanol 0.7744 -> 0.8062
    (fixture) and 0.7963 -> 0.8784 (Crepin's must) against an untouched 0.80; band 0.2256 ->
    0.1938 against an untouched 0.20; **both Minebois legs 1.73x/1.68x -> 1.01x/0.98x her
    published shares**. Uptake draws **no precursor**, so the one channel is the biomass
    denominator. Suite **2015 + 6 xfail**. **D-120's DIRECTION leg is BACK** (D-245's `Flags`
    reverses); its INSTRUMENT leg is only thinner (cap bite 12.7 % -> 4.82 %), so the fifth stays.
  - **D-247 is REINFORCED, not overturned:** both gate rescalings are now **inert** (+1.4e-5,
    +4.3e-6). Do not revisit them.
  - **REFUSED on the way: a storage-quota / NCR feedback.** Crepin's yeast consume 180 mg N/L and
    build ~3 g/L, which at `f_N ~ 0.06` **is** 180 mg N — there is **no extent surplus to bound**;
    the phenomenon is timing. Such a knob is inert where it can be scored and load-bearing only
    where nothing measures it. Do not build it.
  - **RESIDUE, `Flags: D-100`, owner's call:** MLF/Brett draw **only** the amino-acid pair and
    cannot read the `N` slot uptake fills, so the model **over-states yeast/bacteria nitrogen
    competition**. Those tests **isolate the competitor** (the D-33 re-route precedent) — never
    disable the Process globally and never put it behind a scenario flag. `MaillardBrowning`'s
    N-park now reads exactly 0.0 without lees: that EXTENDS D-100/D-104's recorded
    "autolysis-sourced" position from the precursors to the pair; it is not a new decision.

- **`yan_mgl`'s two meanings — MEASURED and REFUSED at D-248. Do NOT re-propose it as open, and
  do NOT ship the compile-seam half alone.** The conflation is **one species wide** on this
  registry (arginine, 4 N vs 3 assimilable; every precursor is 1-N, glutamine's two both release,
  and trp/his are not pools). **The frames CANCEL as an identity**: the seam carves out at TOTAL
  nitrogen and every deamination releases TOTAL, so what the run makes available equals what was
  declared, exactly, for any dose — driven against the DECLARATION (250 -> `X0 + 250/f_N` within
  0.04 % at 0.5 g/L). **Repairing the declaration alone makes the OUTCOME worse: +15.28 mg N/L at
  0.5 g/L (6.1 %) and +30.55 at 1 g/L**, all now realised because D-248's uptake consumes it —
  which is why this could not be measured before that landed, and why the two are ONE repair.
  The complete version is **priced, not built**: `draw_assimilable_nitrogen` books arginine 3-of-4
  and parks the fourth as **excreted urea** (elemental `g N/L`, weight 1.0 on `total_nitrogen`),
  confined to that helper, moving five consumers. Urea = the **ethyl-carbamate precursor**, so a
  beat that builds it gets independent fidelity value.

Measurements: `M:\claud_projects\temp\ferment\d243-wine-nitrogen-audit\` and
`d244-wine-n-evaluation-point\` — `FINDINGS-pre.md`, `probe_band.py`, `probe_ceiling.py`,
`probe_partition_stall.py`, `probe_denovo2.py`, `probe_pins.py`, and the `patch*.py` migrations.
D-248: `d248-nitrogen-uncoupling\` — `measure.py` (the 200x sweep), `probe_residual.py`,
`probe_control.py` (the all-ammonium control), `probe_pins.py`, `probe_phe.py`,
`probe_frame.py` + `probe_cancel.py` (the refusal), and the `fix_*.py` / `rewrite_*.py` edits.
