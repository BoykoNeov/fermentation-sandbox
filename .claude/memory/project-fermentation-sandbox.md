---
name: project-fermentation-sandbox
description: "Fermentation Sandbox status, milestones, repo, and where decisions live"
metadata: 
  node_type: memory
  type: project
  originSessionId: e084eace-c954-47ae-9167-4bbeff335946
---

**Fermentation Sandbox** — research-grade wine/beer fermentation simulation
engine in Python (uv, scipy/numpy/pydantic). Public repo:
https://github.com/BoykoNeov/fermentation-sandbox (default branch `main`).

**Status (as of 2026-06-20): Milestone 0 complete** — the honest, tested skeleton
before any real kinetics. Layered packages (parameters/units → core → runtime →
scenario/validation), tier system, provenance-enforced parameter store, plain
numpy state vector + Process/ProcessSet, SI-ish units (g/L, K, hours), solve_ivp
(BDF) runtime, conservation + benchmark validation harness. 69 tests green; ruff +
mypy(strict on src) clean; CI on Python 3.13/3.14.

**Milestone 1 (in progress):** single-strain, isothermal, nitrogen-limited primary
fermentation that passes the wine (~24 °Brix → dry in 10–14 d) **and** beer
(~1.048 OG → ~1.010 in 5–7 d) §2.2 benchmarks **in parallel**. No real datasets
yet (validate vs published curves; harness is data-ready). Biggest effort =
parameter sourcing.

**M1 progress (2026-06-21):** four tasks done — (1) medium state schemas
(`fermentation.core.media`: `wine_schema` S(1), `beer_schema` S(3), `Medium`,
`MEDIA`, `get_medium`); (2) the scenario→core seam
(`fermentation.scenario.compile_scenario` → `CompiledScenario`; the only
industry→canonical unit boundary); (3) the carbon/nitrogen/mass conservation
quantity functions — new `fermentation.core.chemistry` (molar masses, carbon
fractions, Gay-Lussac split = single source of truth shared with kinetics) +
`fermentation.validation.total_carbon`/`total_nitrogen`/`total_mass`. Biomass C/N
composition is now a Parameter (`biomass_C_fraction`/`biomass_N_fraction` in
`wine_generic.yaml`), passed into the carbon/N builders. New decision **D-8**:
carbon is the rigorous cross-medium invariant; M1 uses the theoretical Gay-Lussac
split (glycerol/realised-yield deferred to Tier-2); mass is wine/hexose-only
(`total_mass` rejects beer's multi-component sugar — hydrolysis water breaks
closure). Then (4) the first kinetic Process, **`GrowthNitrogenLimited`**
(`fermentation.core.kinetics`): Monod growth co-limited by sugar and YAN; draws
N *and* the biomass carbon skeleton from N/S so carbon+nitrogen conserve on its
own to ~1e-13 (no anabolic CO2 in M1 → `total_mass{S,E,CO2}` intentionally open
under growth; biomass N-capped at X0+N0/f_N). `chemistry.sugar_species` added
(S-slot→species map, core can't import validation). Kept OUT of `MEDIA` until the
process set is complete (would break the no-kinetics baseline test). D-8 extended
with the biomass-carbon-routing note. 113 tests green; pushed to `main`.
**M1 batch 2 (2026-06-21):** fifth task done — **`SugarUptakeToEthanolCO2`**
(`fermentation.core.kinetics.uptake`), the fermentative flux. Design forks (user
chose, recorded as **D-9**): uptake is **biomass-catalysed and decoupled from
growth** (`r = q_sugar_max·X·S/(K_sugar_uptake+S)`) — not Pirt growth-coupled —
so it keeps fermenting after N runs out and reaches dryness; beer's slots consume
in preference order via a **smooth catabolite-repression** factor
(`Π_{j<i} K_rep/(K_rep+S_j)`), not a hard switch (kink hurts BDF). Theoretical
Gay-Lussac yields via new `chemistry.HEXOSE_UNITS`/`ethanol_yield`/`co2_yield`
(generalise `*_PER_HEXOSE` to di-/trisaccharides) → carbon (wine+beer) and mass
(wine) close to machine precision; `Y_ethanol_sugar=0.47` stays the unused Tier-2
glycerol hook. Guards mirror growth (clamp `S_i≥0`, zero when `X≤0`). New
speculative params `q_sugar_max`/`K_sugar_uptake`/`K_repression` in
`wine_generic.yaml`. 122 tests green. **Architecture note for next task:**
`ProcessSet` is purely additive, so `EthanolInhibition` can't multiply onto uptake
as a summed Process — it lives inside the rate or in the modifier-hook the
Arrhenius task adds; uptake's rate isn't yet exposed as a wrappable unit, so
expect to touch its structure then (D-9).

**M1 batch 3 (2026-06-21):** sixth task done — **`EthanolInhibition`**
(`fermentation.core.kinetics.inhibition`). Multiplicative, so NOT a summed
`Process`: this task introduced the **`RateModifier`** abstraction
(`fermentation.core.process`). A modifier declares `name`/`tier`/`modifies` (target
Process names)/`reads` and returns a scalar `factor`; `ProcessSet` multiplies that
onto each targeted Process's **whole contribution vector** before summing — which
preserves conservation (uniform scaling of a conserving flux) and needed **no uptake
refactor** (wraps at the set level, not inside the rate — cleaner than D-9 imagined).
Modifiers share the Processes' enable/disable namespace and feed `tier_of` (a
speculative modifier drags its target vars' tier down). Form: Levenspiel/Luong wall
`f = (1 - E/E_max)^n`; `E_max = ethanol_tolerance` (existing param, wall semantics),
new speculative `ethanol_inhibition_exponent` (n=2 → C¹-smooth touchdown, no BDF
kink); targets `SugarUptakeToEthanolCO2` only (growth is N-shut-off before E climbs).
New decision **D-10** (records the deviation: hook built here, one task early;
Arrhenius will reuse it). Known tuning tension: placeholder `E_max=110 < E_final`
(~124-135 g/L), so an inhibited wine run stalls short of dryness — a sourcing/tuning
fix, not a form flaw; tests assert the *mechanism* (smooth / monotone / `[0,1]` /
conservation-preserving / togglable, + a 3-way growth+uptake+inhibition carbon&N
closure), never dryness-under-inhibition. Held OUT of `MEDIA` with the other
kinetics. 142 tests green.

**M1 batch 4 (2026-06-21):** seventh task done — **parameter-tier propagation
(D-1 closed).** `Process` gained `reads: tuple[str, ...] = ()` (now matching
`RateModifier`). `ProcessSet.tier_of`/`tier_map`/`overall_tier` take an optional
`param_tiers: Mapping[str, Tier]`; for each contributing Process *and* each modifier
scaling it, the lowest tier of its `reads` is folded in via `combine`, so a VALIDATED
mechanism on a speculative parameter reports speculative. Two honesty guards: a
declared `read` absent from `param_tiers` **raises `KeyError`** (no silent default to
validated); `param_tiers=None` gives the narrower *structural* (Process/modifier-only)
tier. Threaded end-to-end — new `ParameterSet.tier_map()` (`{name: Tier}`, the tier
counterpart to `resolve()`) → `simulate(..., param_tiers=...)` → `Trajectory.tier_map`.
Advisor steered two calls: (a) the real gap is the runtime path, not just `tier_of` —
a green unit test alone leaves users seeing borrowed credibility; (b) test the actual
`ParameterSet.tier_map()`→`tier_of` bridge, which every other test stubs with a
hand-built dict. D-1 status flipped M0→M1 closed; D-10 note's stale "next task" line
fixed. 149 tests green; ruff + mypy clean.

**M1 batch 5 (2026-06-21):** eighth task done — **`ArrheniusTemperature`**
(`fermentation.core.kinetics.arrhenius`), a `RateModifier` reusing the D-10 hook (no
new mechanism). Built **per rate, not one shared E_a** (advisor: this is the
*established* design, not a disagreement to surface — `E_a_growth` param name +
"A + E_a per rate" context doc already commit it; "targets both" describes the
mechanism, not instance count). Two instances via classmethods:
`ArrheniusTemperature.for_growth()` (reads `E_a_growth`) and `.for_uptake()` (reads
`E_a_uptake`), sharing one `T_ref`. **Reference-anchored** form
`f = exp(-(E_a/R)(1/T - 1/T_ref))`: `f=1` at `T_ref` so the measured rate constant is
used unscaled — **no separate pre-exponential A** (it'd double-book what
`mu_max`/`q_sugar_max` already encode). Always positive (exp) → **no clamp**, and
conservation survives **stacking** (uptake × inhibition × Arrhenius = one scalar on a
conserving vector). Gotchas the advisor flagged: `name` set **per-instance** in
`__init__` (else `ProcessSet` duplicate-name clash); gas constant `R` lives **in the
module** (SI-exact, cited), not `chemistry.py` (stoichiometry-scoped) nor the store
(empirical only); reads `T` from **state** (Kelvin) → already non-isothermal-ready.
New speculative params `E_a_uptake`, `T_ref` (set equal to growth pending sourcing —
separate params are the point; higher E_a = *steeper* T-dependence, so don't infer an
ordering from placeholders). New decision **D-11**. Also closed the **directional
check** task item (warmer-ferments-faster test). 162 tests green; ruff + mypy clean.
Advisor's post-commit catch (fixed): the tier-drag test passed for the wrong reason
(uptake's own reads already speculative) → rewrote it to hold uptake's reads VALIDATED
so only the modifier drives the speculative cap.

**M1 batch 6 (2026-06-21):** ninth task done — **parameter sourcing (D-12).**
`wine_generic.yaml` rewritten from the keystone **Coleman, Fish & Block 2007**
(`10.1128/aem.00670-07`, PDF read directly; strain Premier Cuvée = EC-1118): `mu_max`
0.095/h, `K_n` 0.0088 g/L, `q_sugar_max` **0.85** g/g/h (= β_max/Y_E/S per their eq 5,
*not* β_max — advisor caught this), `K_sugar_uptake` 10.3 g/L, all evaluated at the
20 °C T_ref (Coleman's "Log" is natural log — verified). `E_a_growth/uptake` ≈55 kJ/mol
**derived** from the log-linear T-slope (`E_a = a1·R·T_ref²`; Coleman is polynomial-in-T,
not Arrhenius). `ethanol_tolerance` 142 g/L from the Premier Cuvée TDS (18% v/v) —
resolves the D-10 stall, sourced independently of the benchmark. `beer_generic.yaml`
**created** from open-access **Zamudio Lara et al. 2022** (`10.3390/foods11223602`,
Tables 5/6): `mu_max` 0.098/h, `K_sugar_uptake` 12 g/L; rest transferred/derived.
Tiers promoted to `plausible` **only where a source measures our form** (advisor:
don't inflate) — `K_s`, `K_repression`, `ethanol_inhibition_exponent` stay speculative;
beer is honestly thinner (Droop/Monod-on-sugar models ≠ our Monod-on-N). Beer `E_a`'s
keep the Coleman-derived value, NOT de Andrés-Toro's reported ~35 kJ/mol (primary
paywalled/unread — don't mint a DOI onto an unread number; inert at isothermal M1
anyway). Fixed wrong DOI in context doc (`00845-07`→`00670-07`). Added
`test_load_shipped_beer_parameters`; updated stale "beer file doesn't exist" comments.
163 tests green; ruff + mypy clean.

**M1 batch 7 (2026-06-22): MEDIA wiring + ethanol-brake redesign (D-13).** Wired the
validated-core kinetics into `MEDIA` (`Medium` gained `modifier_factories`;
`build_process_set` threads `modifiers=`). Then the big one — **retired the Luong
ethanol wall from the default media** and replaced it with **`EthanolInactivation`**
(`fermentation.core.kinetics.inactivation`), Coleman's ethanol-driven cell-death term
(eqs 2/7: `r = k'_d·E·X` moves viable `X` → new `X_dead` pool). **Two-pool** biomass
(chosen over a φ-fraction): identical composition ⇒ the X→X_dead transfer is carbon/
nitrogen-neutral *by construction* (`conservation.py` weights `X_dead` like `X`). Why:
the wall is instantaneous+reversible and stalls short of dryness; only a *cumulative*
viability loss sets a wine's timescale and still finishes. **`EthanolInhibition` kept as
an optional class, just unwired.** New decision **D-13**.
**The k'_d archaeology (advisor caught my misframe):** Coleman Table A2 prints the
`Log(k'_d)` `a1` mean as `−1.08e-3`; its 95% CR `[−1.94e-1,−3.30e-2]` is centred on
`−1.13e-1`, and the corrected `−1.08e-1` sits ~dead-centre while the printed `−1.08e-3` is
~1.4 half-widths past the upper bound (opposite side) — a **published typesetting typo**
(dropped `×10ⁿ` exponent). True value `−1.08e-1`. Confirmed 3 ways: reproduces the paper's stated ~13× rise 11→35 °C
(as-printed gives 191×); keeps k'_d(35 °C) under Fig 3b's axis (as-printed overshoots
30×); the IDENTICAL defect sits in the `Log(Y_X/N)` row. Corrected `k'_d(20 °C)=4.28e-5`
(was `3.64e-4`, which stalls even Coleman's *own* reconstructed model at ~108 g/L). PDFs
now in `manual_sources/`; the text layer (decoded `/H110xx` glyphs) settled it. Wired
wine now ferments to dryness. 177 tests green (+10 in `test_kinetics_inactivation.py`);
ruff + mypy clean.

**M1 batch 8 (2026-06-22): N-dependent biomass yield + window re-anchor (D-14).** Task #7
("calibrate the Fig 6c reconstruction") resolved — and **it overturned its own premise.**
The D-13 gap-(a) "uptake-speed gap (~2× too fast)" was a **misdiagnosis**: a faithful
re-implementation of Coleman's comprehensive model (eqs 1-8, Table A2 @ 20 °C — the model
the paper validates against Fig 6c) reproduces our engine **line-for-line** at BOTH 80 and
330 mg N/L (RMSE ~1.3 g/L ≈ 0.5 % of 264 g/L; tracked in `test_coleman_reconstruction.py`).
A 24-Brix/20 °C ample-N wine *should* finish ~6-9 d; the 10-14 d figure was a generic
handoff heuristic, never Coleman. The one **real** gap: a **fixed** `Y_X/N` (= 1/f_N = 8.77)
built too little biomass at low N and **stuck** (~31 g/L residual) where Coleman finishes.
Fix = adopt Coleman's Fig 4 regression `ln Y_X/N = 3.50 − 3.61e-3·YAN_mgL` (the `a1` exponent
carries the **identical published typo as `k'_d`**, D-13: printed −3.61, CR forces −3.61e-3;
reproduces Fig 4 at 80→24.8, 330→10.1 g/g). **Computed at the compile boundary** (new
pattern): unlike the T-regressions pre-baked into YAML at a fixed `T_ref`, this evaluates at
the scenario's *initial* N, so `compile_scenario` overrides `biomass_N_fraction` from the
scenario YAN — keeping the single-source N-conservation contract exact (`total_nitrogen`
reads the same per-run constant; `d/dt[N+f_N·X]=0` for any f_N). **Beer keeps static f_N**
(override gated on the wine-only regression coeffs being present). **Window re-anchored
10-14 → 8-14 d** (user's call on the guarded §2.2 spec); fixture YAN=80 + pitch 0.25 g/L
(both Coleman-anchored, NOT swept to fit — pitch is a real ~2.6-d lever at low N) lands
~9.2 d; benchmark unskipped + passing. Advisor's task-complete pass tightened three things:
Coleman RMSE guard 6.0→2.0 g/L (observed+50 % margin, was 4.5× loose); **asserted the Fig 4
endpoints** (24.8/10.1) directly in `test_compile.py` (the typo correction now has a data
guard, not just prose); inline `Uncertainty` bracket documented as derived metadata; and the
ABV realised-yield/glycerol-sink work promoted from a dangling caveat to its **own tracked
unchecked task** in `milestone-1-tasks.md`. 184 tests green; ruff + mypy clean.

**M1 batch 9 (2026-06-23): beer + CO2 §2.2 benchmarks pass (D-15).** Both unskipped & green.
**Apparent gravity is the crux:** "FG 1.010" is a *hydrometer* (ethanol-depressed) reading,
not real extract (a 1.048 ale's real FG is ~1.016), so model `(sugar,ethanol)` → gravity goes
through the **Balling/Tabarie** relation `RE=0.1808·OE+0.8192·AE`, added to `units/convert.py`
(`real_to_apparent_extract`/`apparent_gravity`, cited). This is load-bearing fidelity: vs
real-extract gravity 1.010 would need ~79 % RDF; apparent lets a realistic all-malt wort hit
1.010. **No new state/param** — unfermentable extract is implicit (`OG_extract − S0`); the wort
spec (S0≈88 g/L — Zamudio's *measured* fermentable sugar — of ~125 g/L total at 1.048, RDF
~70 %; glucose/maltose/maltotriose ≈15/62/23 %; YAN 200, pitch 0.6) lives in the **test
fixture**, sourced NOT back-solved (D-14 discipline); finishes apparent ~1.007 (~3.5-pt margin
below 1.010, robust vs an `inf` flip) → **~5.5 d** (window 5-7). **`q_sugar_max` re-derived
1.5 → 0.5** (the 1.5 = Zamudio `k_S·mu_max` is the growth-COUPLED peak; our decoupled uptake's
realised-equivalent is ~3× lower), still **speculative** `[0.3,1.5]`. **Honesty split (D-15,
advisor-tightened):** the *endpoint* (FG ~1.007, ABV ~5.8 %) genuinely *falls out* of the
sourced wort (real validation); the *timescale* is set by a speculative q at the low end of its
range → confirms q is *consistent with* 5-7 d, NOT that the window emerges unforced — weaker
than wine ("beer is honestly thinner", D-12). **CO2 ratio = 0.977** (not 1.0): the ~2-3 %
deficit is sugar carbon growth routes to biomass (no anabolic CO2 in M1); window `[0.95,1.05]`
accommodates it; test also asserts d(CO2)/dt rises-then-tails.
Added `test_apparent_*` to `test_units.py`. 189 tests green; ruff check + mypy clean.
**Repo-hygiene flag (RESOLVED in commit 81e09bd):** the ruff-0.15.18 format drift across 7
files was swept; `ruff format --check .` is clean again on `main`.

**M1 batch 10 (2026-06-23): realised ethanol yield — glycerol/byproduct sink + must
fermentable fraction (D-16).** Closed the last modelling gap: wine ABV 16.9 % → **~15.0 %**,
realised `Y_E ≈ 0.482`. User chose (over advisor's simpler anchor-to-0.47) **explicit glycerol
+ small byproduct lump AND fix sugar-loading now**; numbers **fall out of independently-sourced
inputs**, not reverse-engineered. (a) Two new produced-only state pools `Gly`/`Byp`, carbon-
accounted as glycerol C₃H₈O₃ / succinic acid C₄H₆O₄. (b) The realised-yield split is **folded
into `SugarUptakeToEthanolCO2`'s yields** — `dS` unchanged (dryness preserved), theoretical
ethanol/CO2 scaled by `(1 − f_C/c(species))`, diverted carbon deposited to Gly/Byp; **carbon
closes for any yields by construction**. Both yields **default 0 = theoretical core** (togglable
per directive 3); **beer carries 0** so its CO2 benchmark is unaffected. (c) `Y_glycerol_sugar`
0.035 g/g (→~8.6 g/L, plausible — Waterhouse Lab/Scanes/Ribéreau-Gayon), `Y_byproduct_sugar`
0.014 g/g (→~3.4 g/L, speculative lump), **`must_fermentable_fraction` 0.93** applied at
compile (only ~93 % of Brix solids are fermentable hexose → loads ~245 not 264 g/L, plausible).
**`VarSpec.default`** added: produced-only pools (X_dead, Gly, Byp) omittable from `pack()` →
2 state vars landed without touching ~37 call sites; substrate vars stay required.
**`total_mass{S,E,CO2}` now asserted only byproduct-off** (glycerol draws redox H like biomass);
`total_carbon` weights Gly/Byp and stays the invariant. **Coleman reconstruction fixed** to feed
the reference the same fermentable S₀ (still RMSE ~1.3 g/L — model-vs-model, shared init).
**Dryness tightened 9.2 → 8.33 d** (still in [8,14]; advisor: *surfaced, NOT tuned* — tuning a
sourced yield to buy margin is what D-14/D-15 forbid). New realism guard
`test_wine_abv_and_glycerol_are_realistic`. Decision **D-16**. 193 tests green; ruff + mypy clean.
Committed `28d3f4e` and **pushed to `main`** (user reaffirmed "always commit and push" — don't
pause to confirm; see [[feedback-batch-end-ritual]]). **Caveat:** ABV 14.98 % sits at the top of
[14,15] — a future Y_E nudge would exit the band.

**M1 batch 11 (2026-06-23): tier-promotion sweep — promote NOTHING (D-17). M1 COMPLETE.**
The final M1 task. Reviewed every Process/modifier/parameter now all three §2.2 benchmarks
pass. **User's call (AskUserQuestion, over a "promote 3 wine mechanisms" alternative): promote
nothing; record why.** Bar: **VALIDATED is reserved for validation against independent *measured*
time-series** (D-C: none exist yet). Passing §2.2 *confirms PLAUSIBLE is earned* (sound forms,
sourced params, reproduces Coleman) but the wine window is re-anchored to Coleman — the same
source the params come from — so it's a **faithful-implementation cross-check, not independent
validation**. **Advisor's linchpin (verified empirically on the real compile path):** promoting
the 3 core Processes to VALIDATED is **near-inert** — *zero* output-tier change on the param-aware
path (D-1 guarantee, what reporting uses) for both media; on the structural path only `X_dead`
flips (it's the lone default-media var no Arrhenius modifier covers). Wine flux vars are param-
capped: `X`/`S`/`N` by `K_s`, `E`/`CO2`/`Gly`/`Byp` by `K_repression`+`Y_byproduct_sugar` (all
speculative). So the tier system already reports honestly regardless of the mechanism-axis label —
parameter-tier propagation (D-1) + Arrhenius capping (D-11) working as designed; the decision is
documentation, not credibility. **Clean calls (hold regardless of bar):** Arrhenius stays plausible
(inert at `T_ref` `f=1`, §2.2 never exercises it — can't promote the untested); beer `q_sugar_max`
+ `K_s`/`K_repression`/`Y_byproduct_sugar`/`ethanol_inhibition_exponent` stay speculative (D-12/15).
Rewrote the pre-registered "promote once §2.2 passes" docstrings (growth/uptake/inactivation) to
say why; clarified Arrhenius; fixed the retired-and-unwired `EthanolInhibition` docstring (advisor
catch — "not yet validated vs §2.2" was misleading for a class §2.2 never runs; also de-staled its
110→142 g/L `ethanol_tolerance` note). Added a VALIDATED-bar note to `ARCHITECTURE.md` + decision
**D-17**. Future promotion trigger = first independent measured time-series into `ReferenceSeries`.
Docs-only, no tier values/logic touched: 193 tests green, ruff + mypy clean. Committed **f24366d**,
pushed to `main` ([[feedback-batch-end-ritual]]).

**Milestone 1 is COMPLETE** — all `milestone-1-tasks.md` items checked off; all three §2.2
benchmarks (wine ~8.33 d / beer ~5.5 d / CO2 0.977) green; the validated core (growth + uptake +
inactivation + Arrhenius, byproduct sink) is built and tier-honest.

**M2 scoping (2026-06-23): Tier-2 scoped — D-18 + the `milestone-2-*.md` trio.** Two opening calls
made by the user via AskUserQuestion, recorded in **DECISIONS D-18**: (1) **pH is a full proton/
charge-balance solver** (track tartaric/malic/lactic/acetic ± carbonic as carbon-accounted state;
per-RHS `Σ charge = 0` root-find for `[H⁺]`; pH a **derived algebraic pure function**, no `dpH/dt`),
chosen over a tracked-pH approximation — discriminator is **compositionality, not accuracy** (tracked-
pH can only do MLF/SO₂ by *scripting*; charge-balance makes them *emerge*). Resolves handoff §7 open
decision #3. (2) **Byproducts/temperature beat is built FIRST**, inverting handoff §6's "pH first" —
because it closes the last skipped benchmark (`test_lower_temperature_is_slower_but_cleaner`),
exercises the **built-but-untested Arrhenius temperature axis** (inert at the isothermal `T_ref`, so
§2.2 never touched it — D-11/D-17), and is the most self-contained T2 physics (esters/fusels are T+N
driven, **no pH dependency**). Advisor sharpened: don't read the user's 3-item list as build order;
collapse derived-vs-integrated (it falls out of richness); fold three couplings into the pH scope —
**evolved-vs-dissolved CO2** (our `CO2` is the cumulative evolved proxy D-15, not the dissolved pool
carbonic needs), **acid carbon vs the D-16 `Byp`=succinic sink** (don't double-count in `total_carbon`),
**pKa(T)**. Stochastic ensemble wrapper = physics-free runtime layer, parallelizable, no scoping gate.
Docs-only (no code), committed + pushed to `main` ([[feedback-batch-end-ritual]]).

**Next: the byproducts beat** — add `esters`/`fusels` produced-only pools (D-16 `VarSpec.default`
pattern) + additive Ehrlich-fusel / ester-synthesis Processes (no `RateModifier` needed); settle the
carbon-accounting sub-decision (recommend route-from-source, D-16 style); source ester/fusel kinetics
(de Andrés-Toro 1998 / Malherbe 2004 + ester lit, mostly `speculative`/directional per §3.5); unskip
`test_lower_temperature_is_slower_but_cleaner`; record as D-19. Then the **pH charge-balance keystone**
(unblocks SO₂ → MLF → mixed cultures). Full design in `docs/plans/milestone-2-{plan,context,tasks}.md`.
No real datasets yet — when the first measured time-series lands in the data-ready harness
(`ReferenceSeries`/`compare_series`), that's the cue to revisit the D-17 sweep and promote validated.

**Don't restate the design here** — it lives in the repo: `docs/DECISIONS.md`
(all design decisions + rationale), `docs/ARCHITECTURE.md`, and the
`docs/plans/milestone-1-*.md` trio (plan/context/tasks). The original brief is
`docs/FERMENTATION_SIM_HANDOFF.md`. See [[feedback-batch-end-ritual]],
[[reference-claude-best-practices]].
