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

- [ ] **MLF-growth — later composition (decision D-23).** Add a growth Process touching `X_mlf`,
      funded from the amino-acid ledger + autolysis. **All prerequisites now landed** (D-33 fusel
      re-route + deamination; D-34 autolysis refill; D-35 event loop for the post-AF pitch), so the
      remaining block is the *consumer* Process itself (composed with a `pitch_mlf` intervention,
      D-36). The AF nitrogen-exhaustion evidence (D-23) is why it cannot be folded into v1.
- [ ] **Mixed cultures / Brett / sour consortium** — resource competition. (After MLF.)
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
