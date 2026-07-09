# Milestone 2 — task checklist

Build order per DECISIONS D-18: **byproducts/temperature first**, then the pH
charge-balance keystone, then SO₂ / MLF / mixed cultures. The stochastic wrapper is
physics-free and can run in parallel. Keep `pytest` / `ruff` / `mypy` green and the
§2.2 trio passing at every step; never weaken a benchmark.

## Active beat — byproducts + temperature axis

- [x] Add `esters` and `fusels` to `wine_schema` and `beer_schema` as **produced-only
      pools** (`VarSpec.default = 0.0`, the D-16 pattern). Both canonical `compile`
      builders updated to list them (D-16 precedent); no fixture churn. Tier-1 suite +
      §2.2 trio + conservation stay green (193 passed, 1 skipped).
- [x] `FuselAlcoholsEhrlich` Process + unit test — additive, produced-only, monotone in
      T, togglable-off. Ehrlich pathway (amino-acid-derived), N-gated on `N/(K_n+N)`;
      non-monotonic-N simplification flagged ⇒ Process tier **speculative**. Embeds its
      own (steeper) Arrhenius factor via the shared `arrhenius_factor` helper.
- [x] `EsterSynthesis` Process + unit test — additive, produced-only, monotone in T
      (warmth favours esters), coupled to the fermentative-flux Monod shape (shares
      `K_sugar_uptake`). Form **plausible**, placeholder rate params cap output at
      speculative. Both wired into `MEDIA` (wine+beer) as a separate `_BYPRODUCT_PROCESSES`
      tuple (isolable); `test_media`/`test_parameters` ripples handled. 211 green.
- [x] **Carbon-accounting sub-decision — option (a)/a1 LANDED (decision D-19).** Ester/
      fusel Processes route their carbon *out of `S`* (booked as ethyl acetate / isoamyl
      alcohol) and `total_carbon` weights the pools, so esters/fusels are real
      carbon-accounted state. `Byp` double-count resolved by carving higher alcohols out
      of `Y_byproduct_sugar` (wine 0.014→0.012; `Byp` = organic acids/polyols only). Draw
      touches only `S` (never `E`/`CO2`). Carbon closes to 1.1e-13; ABV 14.99 %, Y_E 0.482,
      glycerol 8.49, beer CO₂ ratio 0.975 — all §2.2 guards in band. Both carbon sources
      flagged as bookkeeping stand-ins (amino-acid fusel skeleton; ester ethanol moiety
      already in `E`); fusels carry no CO₂ co-product. 213 green. Full write-up in **D-19**.
- [x] Source + reconcile ester/fusel rate + T-sensitivity parameters; replace
      placeholders; tag tiers honestly (all four stay `speculative`, directional only).
      **Hard constraint held:** each `E_a` > `E_a_uptake` (the load-bearing ordering).
      Sourced via de Andrés-Toro 1998 (read in-source through the open CC-BY Pilarski &
      Gerogiorgis 2022, doi:10.3390/pr10112400) for the ester ORDERING, and Mouret 2015
      (doi:10.1016/j.bej.2015.07.017) + Rollero/Mouret 2014 (doi:10.1007/s00253-014-6210-9,
      owner-provided) for the wine reality. `E_a_esters` 75k→80k, others unchanged.
      **Key finding surfaced:** wine ester *synthesis* is weak/non-monotonic in T and the
      wine *liquid* ester fall with T is largely EVAPORATION — a volatilization sink the
      model omits; for wine the warmer⇒more-aroma direction is carried by FUSELS. Full
      record in **DECISIONS → D-19 sourcing step**. 214 green, ruff + mypy clean.
- [x] **OWNER DECISION (2026-06-30): chose option (B)** — build the volatilization /
      gas-stripping sink *first*, then unskip the benchmark honestly (rather than pass the
      combined-total premise as-was). Landed as **decision D-20**.
- [x] **Volatilization / gas-stripping sink — `EsterVolatilization` (decision D-20).**
      New bookkeeping pool `esters_gas` (headspace) + a Process that strips liquid
      `esters` into it on the CO₂-evolution (fermentative-flux) proxy, first-order in
      liquid ester, **carbon-neutral** (`esters`→`esters_gas`, both ethyl acetate; weighted
      in `total_carbon` like evolved CO₂ → closure stays machine-precision). Activation
      energy `E_a_ester_volatil` set **per medium** (the load-bearing split): wine
      `> E_a_esters` ⇒ liquid esters **fall** with T (Rollero 2014 inversion); beer
      `< E_a_esters` ⇒ liquid esters **rise** with T (de Andrés-Toro warm-ale character).
      Esters-only (fusels far less volatile). Schema 11/13→12/14; isolable in
      `_BYPRODUCT_PROCESSES`. New params `k_ester_volatil`/`E_a_ester_volatil` (both media,
      all speculative). Verified empirically @ 14/20/25 °C.
- [x] **Unskipped & passing `test_lower_temperature_is_slower_but_cleaner`** — rewritten
      honest **per medium** (not a combined total, which would hide the wine inversion):
      both media slower-to-dryness + fewer **fusels** when colder; **beer** liquid esters
      fewer when colder; **wine** liquid esters *more* when colder (the D-20 inversion).
      §2.2 trio + carbon conservation stay green. **222 passed**, ruff + format + mypy clean.
- [x] Record outcomes in **DECISIONS D-20**. Remaining for this beat: update
      `milestone-2-plan.md` + ARCHITECTURE if/when the schema-pool count is referenced there.
- [x] **D-21 follow-up — physical Henry's-law stripping + per-medium SOURCED synthesis E_a
      (owner chose the rigorous unified build, prototype-first).** Advisor reconcile: a
      *sourced* Henry's-law stripping is medium-independent (molecule property, Morakul
      2011), so it cannot invert wine without inverting beer — the direction must live in
      **synthesis** `E_a_esters`, which is per-medium sourced. Replaced D-20's fudged
      per-medium `E_a_ester_volatil` with: (i) a shared **physical** partition enthalpy
      `dH_ester_volatil` = 45 000 J/mol (ethyl-acetate Henry's-law, NIST/Sander; Q10 ≈ 1.8)
      + the `E_a_uptake` gas-flow factor; (ii) per-medium **sourced** `E_a_esters` — beer
      **200 000** (de Andrés-Toro steep), wine **55 100 = `E_a_uptake`** (the mapping for
      *flat* integrated production = Mouret's weak wine ester synthesis). Sourced the
      enthalpy from Morakul 2011 + NIST/Sander; prototyped (wine falls 73/61/50, total flat
      114, gas rises; beer rises 22/72/181), then clean impl + **DECISIONS D-21**. No new
      architecture contract (gas flow = `bare_flux · arrhenius(E_a_uptake)`). 222 green,
      ruff + format + mypy clean. The byproducts beat is now physically + sourced-honest.

## Done — carbon-accounting option (a)/a1 (decision D-19)

Landed as variant **a1** (route from sugar + carve `Byp`), exactly the user's call.
Summary (full record in `docs/DECISIONS.md` → D-19):

- **chemistry.py:** added `M_ETHYL_ACETATE` (C4H8O2) / `M_ISOAMYL_OH` (C5H12O) plus
  their `MOLAR_MASS`/`CARBON_ATOMS` entries (single source of truth for draw + check).
- **byproducts.py:** `_draw_carbon_from_sugar` helper splits the carbon draw across
  sugar slots by carbon content (`d[S_i] -= carbon · s_i / Σ s_j c_j`; serves wine's 1
  slot and beer's 3). Both Processes gain `touches=(…, "S")` and draw their species'
  carbon out of `S` — **never `E`/`CO2`**, so `dX`/`dN`/`dE`/`dCO2` stay byte-for-byte.
- **conservation.total_carbon:** weights `esters` (ethyl acetate) and `fusels` (isoamyl
  alcohol). No `Byp` overlap — see carve below.
- **`Byp` double-count resolved:** wine `Y_byproduct_sugar` 0.014 → 0.012 (drops the
  higher-alcohol share, ~0.0017 g/g); `Byp` is now organic-acids/polyols only. Beer
  needed no carve (`Y_byproduct_sugar` = 0).
- **Caveats documented (code + D-19):** fusel carbon is an amino-acid-skeleton
  stand-in (`N` carries no carbon); ester carbon double-represents the ethanol moiety
  already in `E`; fusels carry no CO₂ co-product. All close the ledger exactly;
  none claim metabolic carbon origin.
- **Tier note:** structural-only `tier_of("S")` drops PLAUSIBLE→SPECULATIVE with
  byproducts on; the param-aware tier users see was *already* speculative, so no
  headline change. Not an a1-vs-a2 discriminator.
- **Tests:** byte-for-byte test strengthened to an exact per-RHS carbon-draw balance;
  added a full-ferment `total_carbon` closure test; the integrated drift test reframed
  to the measured trace bound (X/N uncoupled; S/E/CO2 < 0.5 g/L). chemistry tests cover
  the new species.
- **Verified:** carbon closes to 1.1e-13; ABV 14.99 %, Y_E 0.482, glycerol 8.49, Byp
  2.91; beer CO₂ ratio 0.975 — all §2.2 guards in band. 213 passed, 1 skipped (the
  directional benchmark, still owned by the next item).

## Parallel (physics-free) — stochastic ensemble wrapper ✅ DONE (2026-07-01, decision D-24)

- [x] **`runtime/ensemble.py` — `simulate_ensemble` over provenance bands (decision D-24).**
      Takes the full `ParameterSet` (needs the bands), draws `n_members` triangular
      `(low, value, high)` samples (uniform pluggable, zero-width pinned), integrates each via
      `simulate` on a shared `t_eval` grid, returns an `Ensemble`: deterministic **nominal** +
      surviving **members** + each member's sampled param map + derived `tier_map`. Randomness
      lives **only here** (seeded), core stays pure. Aggregates: `median()`, `percentile(q)`,
      `band(name, low=5, high=95)` (outer P5/P95 bracket), `member_trajectory(i)`. Sampling is
      scoped to the **active** Process set's `reads` (no-op draws avoided; `only`/`exclude`
      filters, e.g. `exclude` the pKa set to pin the D-18 pH anchor). Failed members
      (`success=False` *or* a raising RHS) are caught, recorded and counted; past
      `max_failure_fraction` (0.5) the driver raises rather than return a survivorship-biased
      spread. Constraint-group band-overlap checked (yield partition vacuous; `E_a` ordering
      immaterial) — full record in **D-24**.
- [x] Tests (`tests/test_ensemble.py`, 12): determinism (`nominal` == `simulate` on resolved,
      byte-for-byte); seed-reproducibility (same seed identical, different differs); spread
      tracks input band width; `only=[]` degenerate ensemble (members == nominal); **per-member
      conservation** (mass on a toy, carbon on a real wine scenario, using each member's own
      accounting constants); failed-member accounting + survivorship threshold; active-reads
      scoping; tier map unchanged. 274 green, ruff + mypy clean.

## Ensemble follow-ups — attribution, LHS/Sobol, per-member N ✅ DONE (2026-07-01, decision D-25)

- [x] **Per-member nitrogen conservation** (`tests/test_ensemble.py`). Crown-jewel invariant
      extended from carbon to N. Probed first — closes ~1e-12 using each member's own sampled
      `biomass_N_fraction`; biomass is the only N sink (aa-ledger deferred D-23; fusels route C, not N).
- [x] **Spread attribution** (`analysis.attribute_spread`, `tests/test_attribution.py`, 7). SRC²
      variance decomposition post-hoc from `member_params` (no extra runs); shares roll up by
      parameter tier; `1 − R²` reported as the explicit nonlinear/interaction `unexplained` bucket;
      `method="srrc"` rank fallback. Wine: ethanol spread driven by `k_prime_d`/`q_sugar_max`.
- [x] **LHS/Sobol samplers** (`simulate_ensemble(sampler=…)`). `"mc"` default byte-identical;
      `"lhs"`/`"sobol"` via `scipy.stats.qmc` + inverse-CDF, ~8× lower estimator variance at fixed
      budget, center unshifted. Only varying params take a hypercube dim; Sobol needs power-of-two
      `n_members` (raises otherwise). 288 green, ruff + mypy clean. Full record in **D-25**.

## pH / acid charge-balance solver (keystone, D-18) — ✅ DONE (2026-06-30)

- [x] Track tartaric/malic/lactic as carbon-accounted **wine** state (`wine_schema`;
      acetic folded into the `Byp` succinic read; carbonic scoped out — see #1). Beer acid
      system explicitly deferred (phosphate-buffered, no sourced data).
- [x] `Σ charge = 0` root-find for `[H⁺]` solved in pH-space (`brentq`); `pH` and TA as
      **derived pure functions** of state (`core.acidbase`), smooth/fast for BDF, no
      `dpH/dt`. pKa set as provenance-backed parameters (`acidbase.yaml`, all plausible).
- [x] Resolve the three D-18 couplings: **#1 carbonic omitted** (~0.1 % of buffer below
      pH 4 — scoped); **#2 acid carbon vs `Byp`** = *include-by-reading* (read existing
      `Byp` as succinic-equiv, zero new carbon, double-count closed); **#3 pKa(T)** =
      constant pKa (<0.05 units over 10–30 °C — scoped caveat).
- [x] Per-acid concentrations + measured `initial_ph` wired through `compile_scenario`;
      strong cation **back-solved** (inverse anchoring). Series readout in the new
      top-layer `fermentation.analysis`.
- [x] **ACCEPTANCE GATE:** malic→lactic ΔpH ∈ [0.1, 0.3] on a malic-rich must — lands
      **0.225** (`tests/test_acidbase.py` headline). Plus the emergent Byp-driven pH
      down-drift (~0.067) as a second, unscripted demonstration.

## Later beats (stubs, dependency-ordered)

- [x] **SO₂ — free speciation (molecular fraction), readout-only (decision D-22).** Landed
      2026-06-30. Free SO₂ is a wine-only `so2_free` state slot (dosed via `so2_free_mgl`,
      inert), and `acidbase.molecular_so2` solves pH from the organic acids then returns the
      pH-driven **molecular** (antimicrobial) fraction — `1/(1+10^(pH−pKa₁))`, pKa₁ 1.81
      sourced (Usseglio-Tomasset & Bosia 1984), tier plausible. **Readout-only:** SO₂ is NOT
      in the charge balance (the inverse anchoring makes in-balance vs readout give an
      identical t=0 molecular number — see D-22) and is carbon-free, so dosing it leaves pH
      and `total_carbon` byte-for-byte (isolability test). Gates the textbook curve (6/2/0.6 %
      at pH 3.0/3.5/4.0; ~40 mg/L free for 0.8 mg/L molecular at pH 3.5). **No RHS consumer**
      (antimicrobial suppression wires in with MLF/spoilage). **DEFERRED:** the free/**bound**
      (acetaldehyde-binding) split — acetaldehyde is an unbuilt §3.2 byproduct; and SO₂'s
      back-reaction on pH (additive via an `extra_acids` map if a mid-ferment dose event is
      built). 249 green. Full record in **DECISIONS → D-22**.
- [x] **MLF (*Oenococcus oeni*) v1 — conversion-only (decision D-23). LANDED 2026-07-01.** The
      first RHS consumer of `acidbase.molecular_so2` and of pH. Full record + open-knob choices in
      **D-23 → Resolution**. Summary:
      - `core/kinetics/malolactic.py` `MalolacticConversion`: `r = k_mlf·X_mlf·[malate]/(K_mlf+
        [malate])·g_pH·g_EtOH·g_SO₂·γ(T)`; malic (C4) → lactic (C3) + CO₂ (C1) mole-for-mole, so
        **carbon AND mass close on the existing ledger** (no new conservation code). Touches
        `malic`/`lactic`/`CO2`; tier **speculative**.
      - `X_mlf` **dosed-but-inert catalyst slot** on `wine_schema` (`default=0.0`), dosed via
        `mlf_pitch_gpl`; **explicit** (scales the rate, not folded into `k`) so the later growth
        beat is a clean add-a-Process extension.
      - **Open knobs chosen (all speculative):** temperature = **cardinal optimum** (Rosso 1993
        CTMI, `cardinal_temperature_factor`) not Arrhenius; pH gate = logistic; ethanol gate =
        Luong wall (`ethanol_tolerance_mlf` 110 g/L < yeast's 142); molecular-SO₂ gate =
        exponential at the solved pH.
      - **Isolability (2 layers):** value — zero contribution *before* the pH solve when undosed
        (byte-for-byte core, no wasted `brentq`); tier — compile **disables** the Process when
        `mlf_pitch_gpl ≤ 0` so inert `malic`/`lactic` stay VALIDATED (`tier_of` counts enabled,
        not nonzero).
      - **Acceptance (added, not replaced):** new `test_headline_mlf_raises_ph_emergently` uses
        the no-MLF **control difference** `pH_final(dosed)−pH_final(off)` = **0.1813** ∈ [0.1,0.3];
        the algebraic `test_acidbase` headline (0.225) is retained (proves a different thing).
      - **Emergent ethanol "race-or-stall":** a 24-Brix must (~135 g/L EtOH) crosses the 110 g/L
        MLF wall ~day 4, so MLF must complete in that early window or stall — co-inoculation is the
        only viable mode (post-AF doubly blocked: no event loop *and* EtOH > tolerance). 13 new
        tests; 262 green, ruff + mypy clean, §2.2 trio unchanged.
- [x] **Amino-acid ledger — separate yeast/AF beat (decision D-32). LANDED 2026-07-01.** A
      toggleable `default=0` `amino_acids` wine pool (represented as **arginine**) contributing
      to *both* the carbon and nitrogen ledgers, implemented as a **separate isolable Process**
      (`AminoAcidAssimilation`) — a carbon- *and* nitrogen-neutral **swap** (debit the pool,
      refund sugar carbon + ammonium `N`), biomass untouched, so undosed growth + Coleman stay
      byte-for-byte. **Nitrogen-anchored** rate `ρ = ψ·gate(aa)·f_N·base_dx/y_N`; the N-rich
      arginine representative (mass C:N ≈ 1.29 ≪ biomass ≈ 4.3) keeps the carbon refund ≤ 0.30·ψ
      of growth's draw, so it **never creates hexose** — no clamp needed (advisor's carbon-vs-
      nitrogen asymmetry: carbon over-refund = gluconeogenesis, non-physical; N over-refund =
      deamination, deferred with the fusel re-route). **Correctness crux:** the swap anchors to
      growth's *base* rate, so the wine growth Arrhenius (`for_growth` extra target) and the
      carrying-capacity modifier both scale it — landed **fail-first** (guard tests proven failing
      unscaled: `dS = +0.0279` sugar creation at carrying saturation, arrhenius ratio 1.0 vs
      0.445 — then passing scaled). New per-species nitrogen accounting (`NITROGEN_ATOMS`,
      `nitrogen_mass_fraction`); `shared biomass_growth_rate` helper. New speculative
      `amino_acid_assimilation_fraction=0.5` / `K_amino_acids=0.1 g/L` + `amino_acids_gpl`
      scenario key. Dosed = supplementary-YAN feedback (more biomass), isolability undosed-only.
      11 new tests; **406 green** + 5 benchmark, ruff+mypy clean. Full record in **DECISIONS →
      D-32**. **Still-deferred prerequisites for MLF-growth:** (a) the D-19 fusel Ehrlich carbon
      re-route (needs the deamination branch); (b) an autolytic-peptide source to refill the pool
      post-AF (it is empty at the MLF pitch point, D-23).
- [x] **Fusel Ehrlich re-route — the deamination branch (decision D-33). LANDED 2026-07-01.**
      Prerequisite (a) above. A separate wine-only *swap* (`FuselAminoAcidReroute`) re-sources a
      fraction `g = aa/(K_amino_acids+aa)` of Ehrlich fusel carbon off the D-19 sugar stand-in and
      onto the `amino_acids` pool (arginine), **deaminating** the consumed amino acids' nitrogen to
      ammonium `N` — the deamination branch the re-route was deferred on. It **never touches
      `fusels`** (production stays in `FuselAlcoholsEhrlich`; the two share one extracted
      `fusel_production_rate` so the sugar refund matches the draw exactly, via a shared
      `refund_carbon_to_sugar` — the inverse of `draw_carbon_from_sugar` now used by both the swap
      and the re-route). Carbon + nitrogen close by construction; net sugar `−(1−g)·F_c ≤ 0` (never
      creates sugar). **Separate Process was *forced*** by the beer `touches` contract (can't
      declare `amino_acids`/`N` in the both-media producer). **Not modifier-scaled** (it anchors to
      the fusel rate, which carries its own `E_a_fusels` Arrhenius and is scaled by no
      `RateModifier` — contrast the D-32 swap). Caveat: arginine over-releases N (~4× the real
      leucine→isoamyl ratio), a forced single-species-lump consequence. No new params (reuses
      `K_amino_acids`); disabled with the swap when `amino_acids_gpl ≤ 0` (undosed = byte-for-byte
      core). 9 new tests (`test_fusel_reroute.py`); **417 green** + 5 benchmark, ruff+mypy clean.
      Full record in **DECISIONS → D-33**.
- [x] **Yeast autolysis — the autolytic-peptide source (decision D-34). LANDED 2026-07-01.**
      Prerequisite (b) above. `YeastAutolysis` — the **first consumer of `X_dead`** — refills the
      `amino_acids` pool from dead biomass post-AF (*sur lie*), so the pool a later MLF-with-growth
      model draws on is non-empty. First-order in `X_dead`, temperature-accelerated; **nitrogen-
      anchored**: it liberates the dead-cell nitrogen as amino acids (arginine) and routes the
      carbon-rich remainder to a new carbon-only **`debris`** pool (glucan). **Advisor-decided
      fork:** the excess carbon (dead biomass C:N ≈ 4–11 ≫ arginine's ≈ 1.29, so ~86 % of dead-cell
      carbon can't leave as amino acids) goes to `debris` (cell-wall glucan/mannoprotein — the
      *sur lie* lees), **not CO₂** (which would falsely claim autolysis respires the cell and would
      perturb a benchmarked pool). Carbon + nitrogen close *separately*; the excess-carbon split is
      structurally non-negative (no clamp/kink). **Opt-in** (`autolysis_rate_per_h`, the D-30
      carrying-cap pattern — consumes core state, so off by default; undosed = byte-for-byte core).
      Emergent + verified: amino acids rise from empty post-AF, debris outgrows them. New `glucan`
      species; new wine-only `debris` slot (schema 25→26); new params `k_autolysis`/`E_a_autolysis`
      (speculative). 12 new tests (`test_autolysis.py`), incl. an advisor-ordered three-way
      composition test (autolysis feeds while the D-32 swap + D-33 re-route drain — the actual
      MLF-growth-prerequisite config); **429 green** + 5 benchmark, ruff+mypy clean, §2.2 undosed
      trio unchanged. Full record in **DECISIONS → D-34**.
## Event loop (runtime, decision D-35) — the time-driven mechanism

- [x] **Event-loop driver + temperature ramp (decision D-35). LANDED 2026-07-02.** The runtime's
      first *time-driven* mechanism, built verb-agnostic so temperature scheduling and discrete
      interventions share it. `runtime/schedule.py` `simulate_scheduled` segments a run at
      `ScheduledEvent` breakpoints and restarts the pure `simulate` after each state mutation /
      Process-set reconfiguration / parameter update; carries an **external-flow ledger**
      (`final == initial + Σ flows` across the discontinuities — conservation-as-test held) and
      **min-combines the per-segment tier map** so a late-enabled speculative Process drags its
      variables' tier for the whole run. Isolability: `events=()` is byte-for-byte a plain
      `simulate`. **Owner chose the proper temperature ramp now** (not deferred, not a staircase):
      `TemperatureRamp` drives `dT/dt = temperature_ramp_rate` along the piecewise-**linear**
      schedule; the compile boundary segments only at knots where the slope *changes* (collinear
      knots → one segment; flat/single-knot → none, byte-for-byte core) and mints a VALIDATED,
      un-sampled rate parameter only when it ramps. Always-enabled with a `0.0` isothermal default
      (reasoned deviation from the advisor's disable-gate — same guarantees, simpler). BDF
      integrates the constant-slope `T` exactly (verified to 1e-10). Emergent: a warming ramp
      ferments *between* the cold/hot isothermal bounds (Arrhenius reads the true `T(t)`). 20 tests
      (`test_schedule.py` 9 + `test_temperature_ramp.py` 11); 449 green + 5 benchmark.
- [x] **Discrete winemaking interventions (decision D-36). LANDED 2026-07-02.** On the same driver:
      the verb registry at the compile boundary (`add_dap`/`add_so2`/`rack`/`pitch_mlf`), the
      external-flow ledger's payoff (a mid-ferment DAP dose's emergent H₂S-gate response, D-29), and
      reconciling the compile-time MLF disable-gate with a *later* pitch.
- [x] **Ensemble over a scheduled run (decision D-37). LANDED 2026-07-02.** `simulate_ensemble` gained
      an `events` param and routes through `simulate_scheduled` (was wrapping `simulate`); new
      `CompiledScenario.run_ensemble`. Handles per-member Process-set isolation (snapshot/restore),
      schedule-union sampling scope (mid-run-enabled reads), and a per-member external-flow ledger.
      481 green + 5 benchmark.

## Later beats (dependency-ordered)

- [x] **MLF-growth — dynamic `X_mlf` (decision D-38, resolves the D-23 deferral). LANDED 2026-07-02.**
      New `MalolacticGrowth` Process makes the `X_mlf` catalyst dynamic: it grows O. oeni biomass on
      the `amino_acids` pool (D-32/D-34), and since `MalolacticConversion` is linear in `X_mlf`,
      deacidification *accelerates autocatalytically*. **Nitrogen-anchored, carbon shortfall from
      sugar** (mirrors yeast growth / inverts D-34 autolysis; structurally-positive shortfall, no
      clamp) — the higher-fidelity fork over C-anchored/ADI-deamination (owner away; advisor-picked,
      may revisit). `X_mlf` **promoted to real biomass** (now weighted in both ledgers at the biomass
      fractions), superseding the v1 "carbon-free catalyst" claim; the `pitch_mlf` flow now carries
      bacterial C/N (`test_interventions` updated). Gated on **amino acids alone**
      (`amino_acids_gpl>0`, the swap/re-route gate) in its OWN tuple — that alone prevents a tier
      regression on pitched-but-not-aa-dosed D-23/D-31 runs; NOT gated on the pitch, so
      co-inoculation dominance is **emergent** from the ethanol gate (a high-ABV post-AF pitch is
      arrested; a normal-ABV sequential MLF can still grow), matching how conversion is gated. New
      speculative `mu_max_mlf`/`K_aa_mlf`. Scope: growth-only (no bacterial death/decay yet). 15 new
      tests (`test_mlf_growth.py`), incl. the fail-first autocatalysis acceptance (growth on vs the
      Process disabled). **496 green** + 5 benchmark, ruff + mypy clean. Full record in
      **DECISIONS → D-38**.
- [x] **MLF v2 — benign senescence (`MalolacticSenescence`, decision D-41). LANDED 2026-07-06.** Lifts
      the D-39 v1 tradeoff ("no SO₂ ⇒ bacteria never die"): a new pitch-gated Process gives *O. oeni* a
      small always-on **baseline mortality** (`r_sen = k_senescence_mlf·X_mlf·arrhenius(T)`), so a
      pitched, untreated dry wine slowly declines over **weeks-to-months** into the same `X_mlf_dead`
      pool. **Environment-free** (no pH/ethanol/SO₂ term — the D-39 Luong-wall-wipeout crux reused) and
      **Arrhenius, not γ(T)** (warm accelerates, cold preserves). Separate isolable Process ⇒ the D-39
      SO₂ kill stays byte-for-byte; total mortality is now `r_sen + r_death`. New speculative
      `k_senescence_mlf` (5e-4/h, t½ ~58 d, ~100× below the SO₂ kill), reuses `E_a_death_mlf`/`T_ref`;
      no `brentq` (reads no SO₂/pH). C/N-neutral transfer, `X_mlf_dead` a terminal sink (autolysis
      reads only yeast `X_dead`). §2.2 + the 0.1813 deacidification control-difference unmoved. **552
      green** + 5 benchmark, ruff + mypy clean. Full record in **DECISIONS → D-41**. The MLF arc
      (D-23 → D-31 → D-38 → D-39 → D-41) is complete.
- [ ] **Mixed cultures / Brett / sour consortium** — the volatile-phenol spoilage beat (decision
      D-40). Multi-commit arc mirroring MLF (pathway → growth → death → POF+ yeast).
  - [x] **pt1 — Brett phenol pathway, dosed catalyst. LANDED 2026-07-02.** `BrettDecarboxylation`
        (hydroxycinnamics → vinylphenols + CO2) + `BrettVinylphenolReduction` (vinylphenols →
        ethylphenols); Brett carries **both** enzymes so it spoils POF- wine unaided (the canonical
        funk). Gate = `g_SO₂ · γ(T)` only (Brett is acid- + ethanol-tolerant — no MLF pH/ethanol
        walls); warm optimum (32 °C). 5 new wine slots (`hydroxycinnamics`/`vinylphenols`/
        `ethylphenols`/`X_brett`/`X_brett_dead`), `pitch_brett` verb, `_BRETT_GATED_PROCESSES`
        compile gate, `X_brett`(_dead) in `_LEES_SLOTS`. Carbon closes on the existing ledger.
        11 new tests (`test_brett.py`), incl. the emergent Brett-gated headline + SO₂/rack levers +
        post-AF pitch ethanol-tolerance. All params `speculative`. Full record in **DECISIONS → D-40**.
  - [x] **pt2 — `BrettGrowth`, dynamic `X_brett`. LANDED 2026-07-02.** Nitrogen-anchored on
        `amino_acids`; carbon shortfall drawn from **ethanol** (owner fork) so Brett grows in a *dry*
        finished wine → volatile phenols accelerate **autocatalytically**. No sugar Monod; an
        intrinsic logistic carrying-capacity brake `(1 − X_brett/K)` is the only ceiling (Brett has
        no self-arrest, unlike MLF). An early build blew up under **BDF** (`X_brett`→23, aa→−4.5)
        because the `E≤0` guard had no smooth shadow — fixed with the ethanol Monod `E/(K_E_brett+E)`
        (also physically right: Brett grows on ethanol). Regression pinned under BDF + a
        BDF/RK45/LSODA agreement test. Growth gated at the compile seam on pitch **and** aa dose.
        4 new params, 10 new tests. Full record in **DECISIONS → D-40 (pt2)**.
  - [x] **pt3 — `BrettDeath`, the SO₂-driven kill. LANDED 2026-07-02.** Completes the arc's
        control lever: moves viable `X_brett` → `X_brett_dead` on molecular SO₂ (`r = k_death_brett·
        X_brett·(1−g_SO₂)·arrhenius(T)`), so a sulfite addition doesn't just *pause* Brett (the gate's
        `g_SO₂` already does that) — it **kills** it, and the produced-only `ethylphenols` accrual
        halts. Mirrors `MalolacticDeath` (D-39): **SO₂ alone** drives death (Brett's gate has no
        ethanol/pH term to confound, so this is *directly* correct, not a confounder-fix like MLF's);
        **Arrhenius** temperature, not the cardinal γ(T) (so cold *preserves*, matching why Brett
        survives cold cellars). Carbon/nitrogen-neutral transfer (both pools weighted since pt2 — no
        new ledger code). Pitch-gated in `_BRETT_PROCESSES`/`_BRETT_GATED_PROCESSES` (Brett dies
        whether or not it grew, unlike growth); racking already removes it (pt1 `_LEES_SLOTS`). New
        speculative `k_death_brett`/`E_a_death_brett`. 8 new tests (incl. the integration headline —
        SO₂ crashes a *growing* population: `X_brett_dead` accumulates, `X_brett` falls below the
        dose, ethylphenols end below the un-sulfited control). **535 green** + 5 benchmark, ruff +
        mypy clean. Full record in **DECISIONS → D-40 (pt3)**.
  - [x] **pt4 — `YeastPOFDecarboxylation`, POF+ yeast opt-in + emergent reservoir. LANDED 2026-07-06.**
        Closes the Brett arc + the last M2 physics beat. A POF+ (phenolic-off-flavour-positive) primary
        *S. cerevisiae* carries the cinnamate decarboxylase — the *same* reaction as
        `BrettDecarboxylation` (must `hydroxycinnamics` → `vinylphenols` + CO2, carbon-closing 9 = 8 + 1,
        same routing/species) — but **not** the reductase, so during AF it fills the shared reservoir it
        cannot drain: with no Brett the `vinylphenols` **strand** (`ethylphenols` stays exactly 0 +
        VALIDATED), and a later Brett gets a **head start** on the pre-filled reservoir (the emergent
        yeast/Brett coupling, D-26/D-31 parallel). **Forks (owner-decided):** (1) separate opt-in
        Process, not a strain-flag (isolability); (2) pure-enable key **`pof_positive`** (binary strain
        trait; rate stays YAML), **wholly independent of `brett_pitch_gpl`**; (3) carbon routes from
        `hydroxycinnamics` (forced — same reaction). **Rate flux-coupled** (`EsterSynthesis`/
        `AcetolactateExcretion` idiom: catalyst = viable yeast `X` via `fermentative_flux_shape`, NOT
        `X_brett`), so it runs during AF and stops at dryness; **temperature-flat** (no `E_a_pof`, the
        `AcetolactateExcretion` precedent). **Test-design crux (advisor):** stranding is the PRIMARY
        headline (timing-independent control-difference); the head-start is an **early-time /
        time-to-threshold** claim, since conservation forces the *asymptotic* ethylphenols **equal** (same
        total precursor). New speculative `k_pof_decarb` (2.5e-6 mol/(g·h) → ~49 % of the must
        hydroxycinnamic pool converted during AF). 8 new `test_brett.py` tests + `test_media.py`
        `POF_PROCESSES`. **543 green** + 5 benchmark, ruff + mypy clean. Full record in **DECISIONS →
        D-40 (pt4)**.
- **Remaining §3.2 byproducts** — diacetyl (VDK, the lager rest), acetaldehyde
      (early transient peak), H₂S (N/S-deficiency signal). Owner chose to build these
      **one Process per commit, diacetyl first** (decision D-26).
  - [x] **Diacetyl (VDK) — mechanistic 3-pool "diacetyl rest" (decision D-26). LANDED
        2026-07-01.** Owner chose the **C-full** fidelity target (the α-acetolactate
        reservoir is load-bearing — it makes "crash early ⇒ diacetyl rises" and "warm rest
        clears faster" *emerge*; a 2-pool model reproduces neither) and asked for carbon
        **closer to reality** than a return-to-sugar stand-in, so carbon flows through the
        real species `sugar → α-acetolactate → diacetyl + CO₂ → 2,3-butanediol`, every step
        closing on the existing weighted ledger (draw-from-S ⊕ decarb-via-CO₂ ⊕ neutral
        C4→C4 transfer). Three Processes (`core/kinetics/vicinal_diketones.py`), one per
        commit: `AcetolactateExcretion` (flux-linked reservoir, T-flat) → `Acetolactate
        Decarboxylation` (spontaneous, non-yeast-gated, T-critical, `E_a_decarb` high,
        sourced ordering Haukeli & Lie 1978 / Krogerus 2013) → `DiacetylReduction`
        (enzymatic, gated on **viable X**, **no flux term** so it runs during the rest).
        Wired into **both** media (intrinsic yeast metabolism, isolable but always-on like
        esters); params in a new **shared** `vicinal_diketones.yaml` (the load-bearing decarb
        is non-enzymatic ⇒ medium-agnostic). Shared `draw_carbon_from_sugar`/
        `fermentative_flux_shape` promoted to `carbon_routing.py`. **Emergent + verified
        empirically:** beer final diacetyl 0.195/0.040/0.001 mg/L at 10/18/25 °C; wine
        1.011/0.179/0.001 at 14/20/28 °C — warmer monotonically cleaner; warm = peak-then-
        fall; cold strands diacetyl above threshold with an unconverted reservoir.
        `total_carbon` closes to machine precision on a default compiled run; §2.2 trio
        unmoved. All three Processes **speculative** (only the `E_a` ordering sourced).
        30 new tests; **320 green**, ruff + mypy clean. SCOPE: yeast valine-pathway only —
        MLF/citrate diacetyl deferred. Full record in **DECISIONS → D-26**.
  - [x] **Acetaldehyde — transient ethanol-carbon buffer (decision D-27). LANDED 2026-07-01.**
        Produce-then-reabsorb on the *main* ethanol pathway, reusing the D-26 shape (flux-linked
        production + viable-`X`-gated, no-flux reduction) but with **no middle reservoir** — two
        Processes (`core/kinetics/acetaldehyde.py`), one commit. **KEY FORK (advisor-caught,
        owner-decided):** the D-26 preview said "carbon *draw* stand-in", but acetaldehyde's
        product is `E` itself (the uptake Process already does the full lumped sugar→ethanol),
        so a draw-from-sugar would be a *second parallel pathway* = net-new ethanol inflating
        ABV, scaling with pool turnover. Owner chose the **buffer**: production *borrows* a C2
        slice of ethanol (`d[E] -= r·M_eth/M_acet`), reduction returns it — mole-for-mole C2→C2,
        so carbon closes to machine precision touching **neither `S` nor `CO2`**, and the `E`
        endpoint reconverges to the buffer-off core to relative ~1e-8 ⇒ **§2.2 benchmarks
        preserved to far below tolerance** (all 5 unmoved; a ~1e-4 second-order core drift via
        the `E`→viability brake aside). More faithful (acetaldehyde *is* in-transit ethanol
        carbon), not merely benchmark-safe. **Emergent + verified empirically:** wine peak
        37.5 mg/L @ day 2.7/21, beer 38.2 @ day 1.8/14, both reabsorbed to ~0; warmer clears
        faster (55→37→23 mg/L @ 14/20/28 °C). Wired into **both** media (intrinsic, always-on,
        isolable); shared `acetaldehyde.yaml`. Honest caveats pinned: isolability is
        derivative-level (`S`/`CO2`/`N` drift ~1e-4 via the `E`→viability feedback); structural
        `tier_of("E")` drops PLAUSIBLE→SPECULATIVE (param-aware unchanged, the D-26 `CO2`
        parallel). Both Processes speculative. 21 new tests; **342 green**, ruff+mypy clean.
        SCOPE: metabolite only — the SO₂ free/bound split it unlocks is a separate readout beat.
        Full record in **DECISIONS → D-27**.
  - [x] **SO₂ free/bound split — total conserved, free/bound/molecular derived (decision D-28).
        LANDED 2026-07-01.** The readout unlocked by acetaldehyde becoming real state (D-27).
        **Owner-decided fork:** the dosed slot is reinterpreted as **total** SO₂ (rename
        `so2_free`→`so2_total`, conserved/inert) and free/bound are DERIVED at the solved pH by
        the acetaldehyde-bisulfite binding equilibrium (`bound_so2_molar` solves `(A−x)(C−x)β −
        Kx = 0` referenced to bisulfite; `free = total − bound`, `molecular = free × fraction`).
        Option chosen for **conservation** (pinning free would make total grow as acetaldehyde
        rises with no dose — incoherent, and no dip). **Emergent + verified:** dosing 50 mg/L
        total, free SO₂ dips 50→0.9 mg/L at the acetaldehyde peak (day 1.7) then recovers to 50;
        `free+bound==total` to machine precision. At acetaldehyde=0 it collapses to D-22 exactly
        (regression anchor). **Readout-only** (the bound-acetaldehyde-protected-from-ADH RHS
        coupling stays deferred, like D-22 keeping SO₂ out of the charge balance). The one live
        consumer — the MLF antimicrobial gate — now reads the derived free-molecular SO₂ (bound
        SO₂ is not antimicrobial), so a small MLF slip appears during the transient binding
        window (`test_so2_dose_suppresses_mlf_in_a_run` updated, documented). New shared param
        `K_acetaldehyde_so2=1.5e-6 mol/L` (Burroughs & Sparks 1973; plausible; basis pinned to
        bisulfite, ≤5 % vs total-free at wine pH). Caveat: acetaldehyde-only binder ⇒ `bound`
        under-estimates ("total" ≈ "free + acetaldehyde-bound"). 7 new tests; **349 green**,
        ruff + mypy clean. Full record in **DECISIONS → D-28**.
  - [x] **SO₂-bound acetaldehyde protected from ADH — the D-28 free/bound RHS coupling (decision
        D-47). LANDED 2026-07-06.** Lands the coupling D-28 deferred as readout-only: alcohol
        dehydrogenase reduces only *free* acetaldehyde (`acidbase.free_acetaldehyde`), so dosed SO₂
        **locks in** acetaldehyde — a sulfited wine strands a residual (~27 mg/L at a 50 mg/L pitch
        dose, ~0.78 mol per mol SO₂) and its free SO₂ stays depressed. Mechanism + magnitude
        **literature-grounded** (Han et al. 2020; K=1.5e-6 = the literature value; the ~0.76×
        degradation slowdown reproduced at sub-stoichiometric field doses). Owner chose **bake-in,
        default-on** (MLF SO₂-gate precedent); this **intentionally retires the D-22/D-28 "SO₂
        readout-only" invariant** for sulfited runs — but undosed runs are byte-for-byte D-27 (exact
        `so2_total>0` guard), no §2.2 benchmark doses SO₂, carbon still closes, and pH is still not a
        charge actor (SO₂ couples only via charge-free acetaldehyde). Emergent downstream: SO₂ dosed
        *during* AF is now only a **partial** MLF brake (stranded acetaldehyde sequesters the
        antimicrobial pool — "bound SO₂ isn't antimicrobial"), and molecular SO₂ nets *down* over the
        run. CAVEAT (speculative): bound acetaldehyde treated inert-to-ADH ⇒ upper bound on
        persistence. No new params. Rewrote 3 SO₂/MLF tests + added a D-47 `test_acetaldehyde.py`
        section (unsulfited byte-for-byte, rate-throttle, post-AF strands ≪ pitch, carbon closes,
        BDF/RK45/LSODA agreement). **606 green** + 5 benchmark, ruff+mypy clean. Full record in
        **DECISIONS → D-47**.
  - [x] **SO₂-induced acetaldehyde over-production — the transient-peak half (decision D-48).
        LANDED 2026-07-06.** Adds `+k_acet_so2_induced·flux·so2_total` to `AcetaldehydeProduction`
        (glyceropyruvic redox pull, Han 2020), a carbon-exact borrow from E, exact `so2_total>0`
        guard. **The task premise was refuted by data:** D-47 protection *alone* already delivers
        1.3–1.5× the field 0.39 mg/mg end-state slope — the end state is capped by the D-28 binding
        equilibrium, not production — so no additive end-state term fits. Owner chose to scope D-48
        to the transient PEAK. Driver = **total** SO₂ (free collapses to ~0 at the peak ⇒ inert
        there). Magnitude sized by a cross-process constraint: the largest k (4e-3) keeping the
        emergent SO₂ MLF brake in its literature partial-brake regime. **610 green** + 5 benchmark,
        ruff+mypy clean. Full record in **DECISIONS → D-48**.
  - [x] **Excreted overflow-pyruvate pool — the first competing SO₂-binding carbonyl (decision
        D-49). LANDED 2026-07-07.** Part 1 of the D-47/D-28 overshoot fix: the 1.3–1.5× overshoot is
        a real missing mechanism — the model routes 100 % of bound SO₂ onto acetaldehyde, but real
        wine shares it with competing carbonyls (pyruvate, α-KG; Jackowetz & Mira de Orduña 2013).
        Built as an **excreted side pool** (D-19/D-26 idiom), NOT acetaldehyde's on-pathway precursor
        (that rework conflates the intracellular flux intermediate with the extracellular excreted
        residual and would make dosed SO₂ *suppress* acetaldehyde — rejected). `PyruvateExcretion`
        draws C3 from `S`; `PyruvateReassimilation` is **flux-linked (co-metabolic, NOT the no-flux
        ADH gate)**, returning to E+CO₂ (C3→C2+C1), so both terms die at dryness and the pool
        **freezes** at the plateau `k_exc/k_reassim` — a persistent residual pegged to end-of-ferment,
        crash- AND duration-independent (30.0 mg/L at 21 and 40 d). A no-flux viable-X gate drained it
        to ~0 (clean ferment ends with yeast viable) — the mid-build mechanism fix, advisor-confirmed.
        Wine-only; both speculative; carbon closes to machine precision; ABV/CO₂ endpoints preserved
        to rel ~4.4e-5 (≪0.1 %). New `test_keto_acids.py` (19) + 6 existing tests updated for the new
        shared `keto_acids.yaml` + `pyruvate` slot (wine 34→35). **629 green** + 5 benchmark,
        ruff+mypy clean. **Next:** D-50 (α-KG), D-51 (coupled multi-carbonyl SO₂ equilibrium — where
        the slope correction lands). Full record in **DECISIONS → D-49**.
  - [x] **Alpha-ketoglutarate — the third excreted keto-acid SO₂-binder, same structure (decision
        D-50). LANDED 2026-07-07.** `AlphaKetoglutarateExcretion`/`AlphaKetoglutarateReassimilation`
        mirror the D-49 pyruvate pair exactly (flux-linked excretion draws C5 from `S`; flux-linked
        co-metabolic reassimilation freezes a persistent residual at dryness). **The one fork:** the
        reassimilation carbon destination. Advisor-caught insight — pyruvate's `C3→C2(ethanol)+
        C1(CO2)` mole-for-mole split is nearly-isolable *because* it happens to equal the Gay-Lussac
        2:1 fermentation carbon ratio, not because "return to E" is inherently safe; routing α-KG to
        succinate/`Byp` (the "more faithful" α-KG-dehydrogenase reaction) instead would divert
        reassimilation *throughput* (~10–20× the residual) permanently away from ethanol, threatening
        the §2.2 ABV/CO₂ benchmarks — and isn't actually more faithful anyway, since α-KG
        dehydrogenase is repressed under the same anaerobic conditions that make α-KG overflow, and
        the real dominant fate is N-coupled glutamate synthesis (unmodelled). **Decision: mirror
        pyruvate's E+CO2 destination, but fix the ratio** — C5 doesn't divide 1:1 like pyruvate's C3,
        so reassimilation returns carbon at the same 2:1 Gay-Lussac ratio (5/3 mol ethanol + 5/3 mol
        CO2 per mole) rather than mole-for-mole. Residual sized lower than pyruvate's 30 mg/L
        (nominal ~20, **measured 20.0 mg/L exactly** on the acceptance run, pyruvate unchanged at
        30.0). New `total_carbon` weighting term (`alpha_ketoglutarate`, own C5 fraction). Wine-only;
        both speculative; carbon closes to machine precision; combined ABV/CO₂ isolability delta
        **measured** rel ~7.3e-5 (≪0.1 %, roughly double pyruvate-alone's ~4e-5, as expected from two
        detours). **CALIBRATION-PENDING flag for D-51:** both residuals are order-of-magnitude
        estimates the multi-carbonyl equilibrium must re-derive against the field slope, not inherit
        as settled — and D-51 must work in **moles** (SO₂ binds molar concentration): α-KG's higher
        molar mass means 20 mg/L is only ~40% of pyruvate's molar contribution despite being 67% of
        it by mass. 17 new `test_keto_acids.py` tests (36 total) + `test_media.py`
        schema-size/`EXPECTED_PROCESSES` updates. **646 green** + 5 benchmark, ruff+mypy clean.
        **Next:** D-51 (coupled multi-carbonyl SO₂ equilibrium reading acetaldehyde + pyruvate + α-KG
        together). Full record in **DECISIONS → D-50**.
  - [x] **Coupled multi-carbonyl SO₂ equilibrium — the actual D-48 overshoot fix, worked in moles
        (decision D-51). LANDED 2026-07-07.** Generalises D-28's single-carbonyl closed-form
        quadratic to N competing carbonyls sharing one bisulfite pool: `bound_so2_molar` now takes
        `(molar_concentration, Kd)` tuples and solves one shared "reactive bisulfite" root via
        `brentq`, each carbonyl's bound share `Aᵢ·h/(Kᵢ+h)` — verified to reduce EXACTLY to the old
        D-28 quadratic at n=1. Wires acetaldehyde + pyruvate + α-KG together, all in **moles**
        (`M_ACETALDEHYDE`/`M_PYRUVATE`/`M_ALPHA_KETOGLUTARATE`) per the D-50 calibration-pending
        flag. Sourced `K_pyruvate_so2` (5.55e-4 mol/L) and `K_alpha_kg_so2` (1.4e-4 mol/L) from
        Burroughs & Sparks (1973) — the same paper whose acetaldehyde Kd (1.5e-6) matches the
        pre-existing `K_acetaldehyde_so2` exactly, a direct cross-check. **The honest finding (the
        task's other half — re-derive the residual ratio against the field slope, don't inherit
        it): D-51 is a real but PARTIAL fix.** Measured end-state acetaldehyde-vs-SO₂ overshoot at
        the nominal 30/20 mg/L residuals drops from D-48's 1.32/1.44/1.53× (at 50/100/200 mg/L SO₂)
        to **1.15/1.32/1.45×** — competition genuinely narrows the gap, concentrated at low dose
        where the finite keto-acid capacity isn't yet saturated. Pushing both residuals to the top
        of their already-sourced literature uncertainty bands (pyruvate 100 mg/L, α-KG 70 mg/L —
        verified the frozen state actually lands there) narrows it further (0.86/1.10/1.29×) but
        still does NOT close it, especially at 200 mg/L — a structural mismatch (finite-capacity
        Langmuir competitors saturate; the field regression is empirically linear across the tested
        dose range), not a value not yet found. Per the owner's own guardrail ("do not force-fit
        beyond the literature-sourced pool ranges") and advisor concurrence, **shipped the nominal
        D-49/D-50 residuals (30/20 mg/L) unchanged** — tuning one ferment's pool to match a
        237-wine ensemble regression would trade documented provenance for a weaker fitted number.
        One genuine regression-adjacent side effect, fixed honestly: the always-on keto-acid pools
        now ALSO compete for bisulfite in `test_malolactic.py`'s SO₂-dosed MLF run (previously only
        acetaldehyde bound), lowering free SO₂ (~21%→~15% of an 80 mg/L dose) and letting MLF edge
        just past halfway converted (~51%, was ~48%) — test band and docstring updated to the
        measured value, not loosened blindly. **650 green** (646 + 4 new D-51 tests) + 5 benchmark,
        ruff+mypy clean. **Deferred:** closing the remaining ~1.1-1.5× gap needs a different
        structure (not more pool mass) — flagged for a future milestone, not blocking M2. Full
        record in **DECISIONS → D-51**.
  - [x] **H₂S — carbon-free produced pool with an inverse-nitrogen gate (decision D-29).
        LANDED 2026-07-01.** N/S-deficiency signal ("rotten egg"): yeast reduces sulfate faster
        than it can fix the sulfide onto nitrogen skeletons, so production is *de-repressed at
        low N* (inverse of the Ehrlich fusel gate). One flux-linked, temperature-flat Process
        `HydrogenSulfideProduction` filling a new carbon-free `h2s` slot:
        `d(h2s)/dt = k_h2s·X·S/(K_su+S)·K_h2s_n/(K_h2s_n+N)`. **Most isolable beat yet:** touches
        only `h2s` (reads X/S/N, writes none), on no conservation ledger, so disabling it leaves
        every other column's RHS byte-for-byte exact (integrated ~1e-7 = solver-mesh only, not a
        physical coupling); no tier headline (writes a fresh pool nothing reads). New shared
        `hydrogen_sulfide.yaml` with a *separate* `K_h2s_n=0.1 g/L` (YAN scale, NOT the growth
        `K_n=0.0088` — that would be a razor-edge gate) and `k_h2s=2e-6/h`; both speculative.
        **Load-bearing empirical check (before the test):** the cross-must cumulative lever is
        **muted** (80/150/300 mg/L YAN → 0.557/0.542/0.527 mg/L) because N is stripped to ~0 by
        day ~1.3 regardless of dose (no residual-N floor, D-23 gap) — so the anchor is the *gate
        direction*: derivative-level rate(low N)>rate(high N), and integrated the low-YAN must
        makes ~1.8× more H₂S by day 1 despite *less* biomass. SCOPE: produced-only (cumulative
        produced, overstates residual; CO₂-stripping sink deferred — the D-19→D-20 precedent).
        15 new tests; **364 green** + 5 benchmark, ruff+mypy clean.
  - [x] **H₂S CO₂-stripping sink — `HydrogenSulfideVolatilization` (decision D-42). LANDED
        2026-07-06.** Lifts the D-29 produced-only overstatement: a flux-linked, first-order
        Henry's-law sink sweeps the volatile `h2s` into a new carbon-free `h2s_gas` headspace pool
        on the CO₂ stream (`-k_h2s_volatil·flux·f_gas(E_a_uptake)·f_part(dH_h2s_volatil)·h2s`), so
        `h2s` is now the µg/L **residual** reality shows and `h2s + h2s_gas` is cumulative produced.
        The **exact ester D-20/D-21 precedent** but **carbon-free** ⇒ *simpler* (both pools on no
        ledger, transfer neutral by construction — **no `conservation.py` change**; the produced-
        total invariant replaces the ester carbon-closure test). **Flux cancels in the residual**
        (`h2s_ss = k_h2s·gate/(k_h2s_volatil·f_gas·f_part)`) ⇒ residual tracks the gate + T, rises
        as N depletes then freezes at dryness. New speculative `k_h2s_volatil=1.0 L/(g·h)` (→ ~99.7%
        stripped, residual 3.73/2.00/0.91 µg/L @ 14/20/28 °C) + **sourced** `dH_h2s_volatil=17.5
        kJ/mol` (Sander −d ln kH/d(1/T) ≈ 2100 K; exothermic ⇒ +sign, Q10 ≈ 1.3). Honest artifact
        flagged: T-flat production + T-rising stripping ⇒ residual *falls* with a warmer ferment
        (unbenchmarked, directional only). Always-on both media in `_H2S_PROCESSES`; params in the
        shared `hydrogen_sulfide.yaml` (medium-agnostic). New `h2s_gas` slot (wine 32→33, beer
        19→20). 8 new sink tests + flipped run-level assertions (h2s→h2s+h2s_gas across
        `test_hydrogen_sulfide`/`test_carrying_capacity`/`test_interventions`). **561 green** + 5
        benchmark, ruff+mypy clean. Full record in **DECISIONS → D-42**. **The §3.2 aroma beat is
        complete; Milestone 2 physics closes.**
  - [x] **Residual-nitrogen floor — opt-in biomass carrying-capacity cap (decision D-30).
        LANDED 2026-07-01.** Closes the nitrogen gap that muted D-29: growth is the sole N sink
        and strips YAN to ~0 by day ~1.3 regardless of dose. New `BiomassCarryingCapacity`
        RateModifier scales growth's whole contribution by a logistic `clamp(1−X/K, 0, 1)` (K =
        `biomass_carrying_capacity`), so biomass saturates below the N ceiling and dose-dependent
        residual YAN survives — with `dN=−f_N·dX` preserved, carbon+nitrogen still close.
        **OPT-IN / disabled by default:** a residual-N floor is a fundamental DEPARTURE from the
        Coleman anchor (which caps nothing; turning it on in default wine breaks the reconstruction
        RMSE gate at 80 *and* 330 mg/L), so it is wired into wine but the compile seam DISABLES it
        unless a scenario passes `carrying_capacity_gpl`. Disabled ⇒ factor 1 + excluded from
        tiers ⇒ undosed wine byte-for-byte the core (exact 0.0) and growth stays PLAUSIBLE; opt in
        ⇒ structural tier drops to speculative (no param-aware headline — growth already reads
        speculative `K_s`). Emergent: H₂S cross-must lever restored (monotone in dose, span widens
        materially vs the muted core) and dose-dependent residual YAN. New speculative
        `biomass_carrying_capacity=2.5 g/L` (author estimate) in `wine_generic.yaml` + optional
        `carrying_capacity_gpl` scenario key (overrides for sweeps). SCOPE: wine-only (beer
        deferred); MLF-unblock is PROSPECTIVE (MLF v1 has no N gate). 16 new tests; **380 green** +
        5 benchmark, ruff+mypy clean. Full record in **DECISIONS → D-30**.
  - [x] **Default-on N redesign — DECLINED (decision D-43, spike 2026-07-06).** A decision-forcing
        spike + a mass-balance argument proved **default-on residual *assimilable* N is
        Coleman-incompatible regardless of mechanism** (two-pool, cell-quota, satiation): Coleman
        builds biomass by ~day 1.3, pinning external assimilable N to ~0 by then for every dose, so
        no biomass-preserving N model can widen the H₂S lever or leave a late-window residual
        without cutting biomass and breaking the Coleman sugar curve. The deferred note's two
        mechanisms have *opposite* Coleman-compatibility (proline split = Coleman-safe but inert;
        residual-assimilable floor = inherently opt-in). Decision: keep the D-30 opt-in cap as-is,
        no refactor. If the H₂S cross-must lever is ever wanted default-on, re-point the *H₂S gate*
        onto a dose-correlated proxy (an H₂S-model change), not the N model. No source change; the
        negative result is the deliverable. Full record in **DECISIONS → D-43**.
  - [x] **MLF-derived diacetyl — O. oeni citrate co-metabolism + bacterial reduction (decision
        D-31). LANDED 2026-07-01.** The real coupling MLF (D-23) unlocks and the deferred half of
        D-26 (which built yeast valine-pathway diacetyl only). Two new *O. oeni* Processes in
        `malolactic.py`: `MalolacticCitrateMetabolism` co-metabolises a dosed `citrate` must input
        (new C6H8O7 species/slot) into α-acetolactate + CO2 feeding the **shared VDK reservoir** —
        so diacetyl *emerges* from the always-on D-26 decarb + reduction, no new diacetyl kinetics
        — and `OenococcusDiacetylReduction` clears diacetyl on the lees (`X_mlf`-gated). **Why a
        citrate pool (load-bearing):** MLF-diacetyl is a late/post-dryness phenomenon, so its
        carbon can't come from sugar (`draw_carbon_from_sugar` no-ops at `S=0`); citrate is present
        independent of sugar. Lumped carbon-closing stand-in `citrate C6 → acetolactate C5 + CO2 C1`
        (6 = 5 + 1; mass gap = untracked acetate/redox); `k_citrate` low so citrate stays mostly
        unconsumed (the trace diacetyl branch — dominant acetate branch omitted, owned caveat). New
        shared `malolactic_environmental_gate` helper (g_pH·g_EtOH·g_SO2·γ(T)) used by both the
        malate conversion and the citrate branch (coupled to citrate's Monod, NOT malate's r, so
        the post-malate peak survives). Emergent: diacetyl lifts ~2.8× the yeast-only baseline into
        the buttery range (~0.28 mg/L), peaks late (~day 5–6), then reduction clears it; carbon
        closes to machine precision. Isolable: un-pitched/citrate-free = byte-for-byte core, citrate
        VALIDATED→speculative only when dosed. Owner forks: citrate = must input; via shared
        reservoir; add O. oeni reduction now. Deferred: bacterial arrest/death + racking event (so
        "SO₂ locks diacetyl in" / permanent stranding), citrate in the pH balance. 14 new tests;
        **395 green** + 5 benchmark, ruff+mypy clean. Full record in **DECISIONS → D-31**.
  - [x] **MLF v2 refinement — bounded ethanol/starvation stress multiplier on `MalolacticSenescence`
        (decision D-52). LANDED 2026-07-07.** With M2 physics all complete through D-51, closes one
        of the two remaining MLF v2 deferrals (`BrettSenescence` stays open — see below). Owner asked
        for "whichever is closer to reality"; an advisor call reversed the first-pass pick
        (`BrettSenescence`) after verifying against DECISIONS D-40 pt3 that Brett's *persistence* is
        already an intentional fidelity choice, not a gap — so a senescence twin would be a downgrade,
        not a gain. `MalolacticSenescence`'s rate gains a bounded stress multiplier:
        `stress = 1 + k_senescence_ethanol_scale·[E/(E+ethanol_tolerance_mlf)] +
        k_senescence_starvation_scale·[K_aa_mlf/(K_aa_mlf+amino_acids)]`, both terms smooth Monod-type
        factors in **[0,1)** (not the Luong wall's near-binary shape that caused the D-39 wipeout),
        capped at 1+1.0+0.5=2.5× **at T_ref** — worst-case half-life ~23 d, verified empirically at the
        RHS level, nowhere near the ~1-week wipeout regime. **A second advisor pass** (post-commit)
        caught that the wipeout test implicitly fixed T=20 °C while claiming "worst case" — split into
        a T_ref-scoped test plus a warm-temperature test proving the invariant that actually holds at
        any temperature (chronic senescence stays far below the acute SO₂ kill, temperature-invariant
        ratio by construction). Reuses `ethanol_tolerance_mlf`/`K_aa_mlf` (no new concentration
        scales); only two new dimensionless ceiling params. Reads no SO₂/pH — still cheaper than the
        SO₂ kill. Genuine side effect measured and banded honestly (the D-51 discipline):
        `test_headline_citrate_lifts_and_then_clears_diacetyl`'s final/peak ratio rose to ~0.861 (less
        viable bacterial reductase late in a 30-d run) — band widened to 0.90 with the measured value
        recorded, not loosened blindly. **OWNER-FLAGGED OPEN QUESTION — RESOLVED in D-53 (see below):**
        `k_senescence_mlf` (5e-4, unchanged at ship time) was D-41-calibrated so a *typical* wine
        loses ~half its O. oeni over ~2 months — but stress≈2× is close to the *typical* case, so the
        typical-wine half-life D-52 actually produced was ~29 d. Owner asked for research rather than
        picking a number; D-53 found real-wine CFU data does not support ANY weeks-to-months
        spontaneous decline, corrected the magnitude ~50x down, not just re-anchored. 4 net new/split
        tests in `test_malolactic.py` + 2 tests re-measured (not weakened) + the diacetyl band update.
        **654 green** + 5 benchmark, ruff+mypy clean.
        **Deferred:** a `BrettSenescence` twin for the D-40 arc remains open, but is
        now framed as a *declined-by-default* option (Brett's persistence is the honest model) rather
        than a straightforward extension — revisit only if a source is found that Brett *does* decay
        benignly in practice. Full record in **DECISIONS → D-52**.
  - [x] **Correction: `k_senescence_mlf` magnitude was wrong by ~50× (decision D-53). LANDED
        2026-07-07.** Owner asked for deep research rather than picking a number to close D-52's
        open calibration question. Finding: real, unsulfited finished wine shows **no detectable
        spontaneous O. oeni decline for 3–5 months** (Windholtz et al. 2025, OENO One, doi:10.20870/
        oeno-one.2025.59.3.9346; Millet 2001 thesis) — the steep decline D-41's citations implied is
        actually **SO₂-driven** (Kioroglou et al. 2020, doi:10.3389/fmicb.2020.562560), i.e. already
        `MalolacticDeath`'s (D-39) territory, not spontaneous senescence. D-41's citations supported
        general "SO₂ controls spoilage LAB" practice, not a specific weeks-to-months claim — a
        misread that propagated uncaught into D-52. Fix: `k_senescence_mlf` 5e-4 → **1e-5** (a round,
        upper-bound-consistent value, not a fit — no source measures decline past 5 months); worst-
        case D-52 combined stress now gives a multi-year half-life. **D-52's stress-multiplier
        mechanism is completely unchanged** — only the baseline it scales was wrong. **Honest
        consequence surfaced to the owner, not buried:** at this magnitude D-52 is now empirically
        inert on every timescale the model simulates. Owner chose to keep the structure (least
        churn) over stripping it back to D-41's flat form. **Test consequence was an assertion
        flip, not a re-band:** `test_so2_crashes_bacteria_over_the_slow_senescence_baseline`
        asserted a *measurable* decline that the new evidence contradicts — renamed and flipped to
        assert near-stability (measured ~0.990 at day 21, was ~0.608); the diacetyl clearing test's
        rationale was corrected too (ratio reverts to ~0.742, closer to D-41's original picture).
        **654 green** (same count, tests reassigned not added), ruff+mypy clean. **Method beat:**
        the third `advisor()` call in the D-52/D-53 arc — each pass caught something the previous
        one missed (wrong pick → test-scope gap → now a magnitude correction the owner's research
        request itself surfaced). Full record in **DECISIONS → D-53**.

- [x] **§3.3 Hop bittering → IBU (decision D-64). LANDED 2026-07-10.** The "additives with clear
      mechanisms" beat, owner-selected off the post-D-63 menu. Two-regime physics in the two places
      each belongs: (1) the **boil** isomerization `alpha --k1--> iso-alpha --k2--> degradation` is a
      wort-side CLOSED-FORM compile-seam calc (Malowicki & Shellhammer 2005, sourced Arrhenius k1/k2;
      doi:10.1021/jf0481296), wired into a new beer-only `iso_alpha` state at t=0 like `initial_ph` —
      NOT a boil ODE phase; (2) the **fermentation** loss `IsoAlphaAcidLoss` (X-gated first-order
      adsorption, the only Process) makes finished IBU fall below the end-of-boil value. **Off the
      carbon ledger** (exogenous hop mass, like dosed SO2) ⇒ `total_carbon` byte-for-byte unchanged
      (asserted directly). **A `hop_utilization_efficiency` (0.55, spec) is ADDED, not fitted** — set
      from literature-typical utilization to avoid the ~2× finished-IBU overprediction of raw kettle
      kinetics (a fidelity failure), with Tinseth as an independent fit-vs-fit *check* (canonical
      recipe ~17.0 IBU vs Tinseth ~17.3), preserving the firewall. New scenario surface:
      `HopAddition` + `hops`/`batch_volume_liters` (the genuinely-new volume quantity)/`boil_celsius`.
      **Tier derives, not asserted:** sourced boil kinetics plausible, capped speculative at the
      finished pool by the loss/efficiency inputs (D-1). Isolable/beer-only: hops on wine is a loud
      error; unhopped beer disables the loss (VALIDATED slot), byte-for-byte the prior core. New
      `hops.yaml` + `core/kinetics/hops.py` + `analysis.ibu_series`; `tests/test_hops.py` (20 tests);
      **704 green**, ruff+mypy clean. **DEFERRED (v1):** gravity-dependence of utilization, dry-hop/
      whirlpool bitterness, hop-form/pH effects, oxidized-alpha (humulinone) bitterness. One advisor
      call before writing (shape endorsed; 5 sharpening points applied). Full record in **DECISIONS
      → D-64**.

- [x] **§3.3 Acid/sugar adjustments (decision D-65). LANDED 2026-07-10 — §3.3 now COMPLETE.** The
      last of the four "additives with clear mechanisms," owner-selected as the natural continuation
      of D-64. Two new intervention verbs at the compile→core seam, **no new Processes** (the brief's
      "simple state mutations via events"), riding the D-35/D-36 external-flow ledger: (1) **`add_acid
      {acid, gpl}`** — general over the D-18 charge-active acids (tartaric/malic/lactic), wine-only by
      slot presence; the pure acid lands on its slot but NOT on `cation_charge`, so the charge balance
      re-solves and **pH drops / TA rises emergently** (potassium-bitartrate deacidification deferred).
      (2) **`add_sugar {sugar_gpl}`** — chaptalization: doses SUCROSE, inverted AT THE DOSE to
      hexose-equivalent (×~1.0526, the hydrolysis-water mass gain) onto wine's `S` or beer's GLUCOSE
      slot specifically (fructose lumped as glucose-equiv — isomers, exact carbon); more sugar ⇒
      higher finished ethanol, emergent. Both book **positive** carbon external flows (the copper
      mercaptan −C mirror), nitrogen-free, ledger closes to machine precision. **Three owner forks
      decided by AskUserQuestion up front** (all chose the more capable option): general acid verb,
      sucrose-with-inversion, wine+beer scope. One new param — `sucrose_inversion_mass_ratio` in
      `additions.yaml`, VALIDATED exact stoichiometry, zero-width band (the `dap_nitrogen_fraction`
      precedent). **No tier movement** (neither verb enables a Process; inert slots). `tests/
      test_interventions.py` +13 (50 total); **717 green**, ruff+mypy clean; every prior benchmark
      unchanged. **DEFERRED (v1):** K-tartrate deacidification, kinetic sucrose/invertase pool, direct
      glucose/fructose dosing, volume change. Full record in **DECISIONS → D-65**.
