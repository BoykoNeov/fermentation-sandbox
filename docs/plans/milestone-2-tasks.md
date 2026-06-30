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

## Parallel (physics-free) — stochastic ensemble wrapper

- [ ] Runtime wrapper over `simulate` sampling each parameter within its provenance
      `Uncertainty` band; run ensembles; aggregate (median + spread). Lives in `runtime`,
      outside the pure core (handoff §1.6). Tier/uncertainty bands reported on outputs.
- [ ] Tests: determinism preserved for a single (unsampled) run; ensemble spread tracks
      the input uncertainty; core stays reproducible.

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

- [ ] **SO₂** — free/bound/molecular equilibrium, pH-driven speciation (pKa ≈ 1.81);
      acetaldehyde binding. (Needs pH.)
- [ ] **MLF (*Oenococcus oeni*)** — second-organism Process on a "pitch MLF" event;
      malic→lactic+CO₂; pH/ethanol/SO₂/T-sensitive growth. (Needs pH.)
- [ ] **Mixed cultures / Brett / sour consortium** — resource competition. (After MLF.)
- [ ] **Remaining §3.2 byproducts** — diacetyl (VDK, the lager rest), acetaldehyde
      (early transient peak), H₂S (N/S-deficiency signal).
