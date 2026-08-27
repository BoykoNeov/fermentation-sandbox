---
name: banded-undrawn-census
description: "D-240/D-241 — the banded AND compile-read AND never-drawn census: 28 names, six classes, priced; SIX are now REPAIRED and drawn"
metadata:
  node_type: memory
  type: project
---

**Live prohibitions — the banded-and-never-drawn census (D-240).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about "parameters the ensemble cannot draw", a seed's uncertainty being
missing from a reported band, `burst_antioxidant_initial`, beer's `*_typical_wort` levels, or
wiring a seed into a Process's `reads`. Every bullet is *what it forbids* + the record to read for
*why*. **If a prohibition looks unconvincing, go read D-240 — do not argue past it from here.**

- **The census is RUN and CLASSIFIED — never re-propose it as unenumerated.** D-237 §6's parked
  item is CLOSED. **32** compile-read names drawn nowhere, **27** banded union-wide / **28**
  per scenario in **six** classes, each with a verdict in `tests/test_banded_undrawn_census.py`'s
  `VERDICTS`. It is
  the **complement** of D-234's census within the compile reads (the two partition them, per
  scenario) and it **OVERLAPS** D-153…D-159's drawability surface — 11 of the 28 are D-159's
  structural class by the other predicate. **Never call the two disjoint**; D-234's disjointness
  assertion is about *its* set and this one, not about D-159's.
- **`copper_typical` is in BOTH registries and that is CORRECT.** Membership is
  scenario-conditional: a D-234 member in an aged wine, a member here in a wine that never ages,
  because `phenolic_browning` is its only reader and is disabled until `begin_aging`. Score any
  such registry **per scenario**; a test that forbade the overlap would go red on the first
  honest one.
- **NEVER quote the UNPAIRED widening (1.672×) — it is noise.** Two ensembles over different
  `only=` sets draw different members. Paired on identical draws the wine ferment's ethanol
  figure is **0.323** [[feedback-pair-the-arms-before-comparing-spreads]]. **The `as shipped`
  arm is EXACTLY 0.000000 on every slot** — forcing all 28 into `only=` moves nothing, because
  none is read at runtime and `y0` was compiled before the draw. That zero is the negative
  control every ratio is read against; do not drop it.
- **The hidden spread is NOT uniform and the 50× band is not the top of it.** Paired hidden/full:
  `methanethiol` **1.88**, `dms` **1.09**, aged `A420` under `direct_burst` **0.705**, ferment
  `E` **0.323**, beer's `X`/`E` **0.000**. Edge-to-edge (upper bound, one name at a time):
  `must_fermentable_fraction` **0.574** of the ethanol spread from a **1.06×** band, against the
  drawn control `mu_max` at **0.441**. **Band WIDTH does not order this table.**
- **Beer's eight `*_typical_wort` are ABSORBED at t=0 to ≤1.2e-13 — never file them as live
  seeds.** The cation anchor is back-solved *through* them (D-238 §4 / D-239 §5 signature). What
  the band buys is a day-1 wobble that decays: acetic ±1.11e-03 → ±9.8e-05 by day 7, formic
  ±9.9e-04 → ±1.4e-06. D-239 already spent that day's envelope on a term with a mechanism.
- **`o2_wort_aeration_beer` is banded 1.45×, live, and worth EXACTLY ZERO** (final `X` moves
  7.7e-09, `E` 5.1e-11). D-213's aging-gate scope decision, re-measured on the trajectory.
  **Never file it as a gap** — it is the zero the other prices are read against, and since D-241
  it is DRAWN and re-seeded per member as the repair's null control: beer's whole reported band
  is unchanged to **1.000** while the rule demonstrably fires (~1e-8 relative).
- **`copper_h2s_binding` / `copper_mercaptan_binding` lose their own `min()` at any real dose.**
  The verb removes `min(pool, copper_gpl · binding)` and a 0.5 mg/L dose is in **126×–2200×**
  excess of the sulfide present, so the constant never enters the arithmetic. **This is D-159's
  supply-limited warning with the sign REVERSED** (excess, not absence) — so the guard pins the
  **capacity ratio**, never the zero: a zero-only test survives a later change to the dose
  arithmetic. `copper_fining_residual_fraction` DOES move (0.068 of the `copper` spread).
- **SIX of the eight live seeds are REPAIRED and DRAWN since D-241 — do not re-propose them.**
  `dms_potential_initial`, `bound_h2s_initial`, `bound_methanethiol_initial`,
  `must_fermentable_fraction`, `o2_wort_aeration_beer`, and `burst_antioxidant_initial` (under
  `direct_burst` only). They MOVED to `tests/test_compile_sampled_census.py`'s `CENSUS`; the
  mirror pin is `test_the_six_repaired_seeds_really_are_drawn_now`. → `seed-reads-repair.md`
- **D-240's decline was an argument against ONE route, not against the repair.** It said a
  repair needs a `y0_for_member` rule **AND** a declaration, and that the declaration is a tier
  claim (`reads` has two masters, D-160). True of declaring on a **Process**. `seed_reads` takes
  the sampling half alone [[feedback-a-blocker-may-be-one-doors-lock]]. **The tier half is still
  unmeasured and is still owed** — never quote D-241 as having settled it.
- **`test_the_priced_names_are_still_undrawn` FORBADE NOTHING against that route** — it scored
  the resolver without `seed_reads`, so it stayed GREEN through the whole repair. D-240 Arm C's
  lesson one layer up: the INSTRUMENT, not just the assertion, must reach the subject
  [[feedback-a-guard-must-be-scored-where-its-subject-lives]]. The shared `_census` passes it now.
- **The two Coleman coefficients are SUBSUMED, not a second beat.** `biomass_N_yield_log_intercept`
  /`_slope` seed no slot — they derive the `biomass_N_fraction` override, and **that parameter is
  itself sampled**, over `[0.03, 0.15]`, which strictly CONTAINS the `[0.051432, 0.108338]` the two
  coefficients imply and is 2.11× wider. Drawing them would double-count one quantity. **Never
  propose a `values_for_member` hook for them.**
- **Do NOT widen the battery to closure/toast/spirit variants.** ~40 more names, zero findings —
  `oak_yield_vanillin_heavy` carries `_medium`'s verdict for `_medium`'s reason. The skip is
  stated, not hidden [[feedback-count-and-print-your-skips]].
- **A registry row must be falsifiable in the direction it is NOT about.** The first
  `test_no_classified_name_has_gone_stale` subtracted `_ZERO_WIDTH` from its own union — the five
  zero-width rows could then never go stale, so the seam could stop reading one and nothing would
  notice. Score them against the **unfiltered** undrawn set instead; verified by planting an
  unreached name in `_ZERO_WIDTH` and watching both tests go RED
  [[feedback-a-guard-must-be-scored-where-its-subject-lives]].
