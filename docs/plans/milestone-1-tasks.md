# Milestone 1 — task checklist

- [x] Define wine `StateSchema` ({X, S(1), E, N, T, CO2}) and beer schema
      ({X, S(3), E, N, T, CO2}). → `fermentation.core.media` (`wine_schema`,
      `beer_schema`, `Medium`, `MEDIA`, `get_medium`).
- [x] `compile(Scenario) -> (y0, ProcessSet, params)` — the scenario→core seam,
      converting industry units to canonical at the boundary. →
      `fermentation.scenario.compile_scenario` returning `CompiledScenario`
      (`y0`, `process_set`, `parameters`/`param_values`, `schema`, `t_span_h`).
      Process set is empty until kinetics land; tests in `tests/test_compile.py`,
      `tests/test_media.py`.
- [x] Carbon + mass conservation quantity functions for the real chemistry.
      → `fermentation.validation.total_carbon` / `total_nitrogen` / `total_mass`,
      weighting state vars by the shared stoichiometry in
      `fermentation.core.chemistry` (single source of truth with the kinetics).
      Carbon + nitrogen are the rigorous atom balances (biomass C/N fraction is a
      passed-in Parameter — `biomass_C_fraction`/`biomass_N_fraction` in the
      store); mass is scoped to the abiotic `S+E+CO2` conversion. Decision D-8;
      tests in `tests/test_chemistry.py`, `tests/test_validation.py`.
- [x] `GrowthNitrogenLimited` Process + unit test. → `fermentation.core.kinetics`
      (`GrowthNitrogenLimited`). Monod growth co-limited by sugar and YAN
      (`mu = mu_max·S/(K_s+S)·N/(K_n+N)`, `dX/dt = mu·X`); draws nitrogen *and* the
      biomass carbon skeleton from `N`/`S` so it conserves `total_nitrogen` and
      `total_carbon` on its own to machine precision (no anabolic CO2 in M1 —
      D-8). Held **out of the `MEDIA` registry** until the process set is complete
      (a growth-only medium is a half-model and would break the no-kinetics
      baseline test). Tests in `tests/test_kinetics_growth.py`. The S-slot→species
      map moved core-ward to `fermentation.core.chemistry.sugar_species` (core may
      not import validation), shared by the Process and the carbon check.
- [x] `SugarUptakeToEthanolCO2` Process (Gay-Lussac yield; vector-aware for beer)
      + unit test. → `fermentation.core.kinetics.uptake`. Biomass-catalysed
      (`r = q_sugar_max·X·S/(K_sugar_uptake+S)`), decoupled from growth so it
      ferments to dryness after N runs out; beer's slots consumed in preference
      order via smooth catabolite repression (`Π_{j<i} K_rep/(K_rep+S_j)`).
      Theoretical Gay-Lussac yields from `chemistry.ethanol_yield`/`co2_yield`
      (new `HEXOSE_UNITS`) so carbon (wine+beer) and mass (wine) close exactly —
      `Y_ethanol_sugar` stays the unused Tier-2 hook (D-8). New speculative params
      `q_sugar_max`/`K_sugar_uptake`/`K_repression`. Decision **D-9**; tests in
      `tests/test_kinetics_uptake.py`.
- [x] `EthanolInhibition` Process + unit test. → `fermentation.core.kinetics.inhibition`.
      Multiplicative, so **not** a summed Process: introduced the `RateModifier` hook
      (`fermentation.core.process`) — `ProcessSet` scales each targeted Process's whole
      contribution vector by the modifier's `factor`, preserving conservation and
      needing no uptake refactor. Levenspiel/Luong wall form `f = (1 - E/E_max)^n`
      (`E_max = ethanol_tolerance`; new speculative `ethanol_inhibition_exponent`, `n=2`
      for a C¹-smooth touchdown), targeting `SugarUptakeToEthanolCO2` only (growth is
      N-shut-off before E climbs). Modifier tiers wired into `tier_of`. Decision
      **D-10**; tests in `tests/test_kinetics_inhibition.py` + modifier-mechanism tests
      in `tests/test_process.py`. Held out of `MEDIA` with the other kinetics; the
      placeholder `E_max=110 < E_final` stall is a tuning-task item, benchmark stays
      skipped.
- [x] **Wire parameter-tier propagation into `ProcessSet.tier_of`** — each Process
      declares the params it reads; an output's tier = lowest of (its Processes'
      tiers, those params' tiers). Closes the D-1 gap so a VALIDATED process on
      speculative params reports speculative. → `Process` gained a `reads` attribute
      (matching `RateModifier`); `tier_of`/`tier_map`/`overall_tier` take an optional
      `param_tiers` map and fold in the lowest tier of each contributor's `reads`
      (process *and* modifier), raising `KeyError` on a declared read missing from the
      map (no silent default to validated). Threaded end-to-end:
      `simulate(..., param_tiers=...)` → `Trajectory.tier_map`, built via new
      `ParameterSet.tier_map()`. Decision **D-1** status flipped to closed; tests in
      `tests/test_process.py` (parameter-tier propagation) +
      `tests/test_integrate.py::test_trajectory_tier_map_caps_on_param_tiers`.
- [x] `ArrheniusTemperature` modifier (isothermal runs for M1) + unit test. →
      `fermentation.core.kinetics.arrhenius`. Reuses the D-10 `RateModifier` hook (no new
      mechanism). **Parameterised per rate** (growth and fermentation differ in T-
      sensitivity): `ArrheniusTemperature.for_growth()` (reads `E_a_growth`, scales
      growth) and `.for_uptake()` (reads `E_a_uptake`, scales uptake) share one `T_ref`.
      Reference-anchored form `f = exp(-(E_a/R)(1/T - 1/T_ref))` — `f=1` at `T_ref`, >1
      above (faster), <1 below; always positive so **no clamp** and conservation is free
      even when stacked with inhibition on uptake. `R` lives in-module (SI-exact, with
      citation); reads `T` from state (Kelvin) so it is T2-ready. New speculative params
      `E_a_uptake`, `T_ref`. Decision **D-11**; tests in `tests/test_kinetics_arrhenius.py`
      (incl. a 4-modifier full-run carbon/nitrogen closure and a warmer-ferments-faster
      directional check). Held out of `MEDIA` with the other kinetics.
- [x] Source + reconcile parameters; replace `wine_generic.yaml` placeholders and
      add `beer_generic.yaml`; re-tag tiers. (Biggest effort — handoff §2.3.) →
      Wine sourced from **Coleman, Fish & Block 2007** (`10.1128/aem.00670-07`, PDF
      read directly): `mu_max` 0.095/h, `K_n` 0.0088 g/L, `q_sugar_max` 0.85 g/g/h
      (= β_max/Y_E/S, eq 5), `K_sugar_uptake` 10.3 g/L, all at the 20 °C T_ref;
      `E_a_growth/uptake` ≈55 kJ/mol derived from the log-linear T-slope
      (`E_a = a1·R·T_ref²`); `ethanol_tolerance` 142 g/L from the Premier Cuvée /
      EC-1118 TDS (18% v/v). Beer added from **Zamudio Lara et al. 2022**
      (`10.3390/foods11223602`, open-access, Tables 5/6): `mu_max` 0.098/h,
      `K_sugar_uptake` 12 g/L; the rest transferred/derived and kept speculative.
      Tiers promoted only where a source measures *our* form; `K_s`,
      `K_repression`, `ethanol_inhibition_exponent` stay speculative. Tests updated
      (`test_compile.py`, `test_parameters.py`). Decision **D-12**.
- [x] **Wire the validated-core kinetics into the `MEDIA` registry.** `Medium` gained
      a `modifier_factories` tuple alongside `process_factories`, and
      `build_process_set` now threads `modifiers=` into the `ProcessSet`. Both wine and
      beer wire the same mechanism set — `GrowthNitrogenLimited` +
      `SugarUptakeToEthanolCO2`, scaled by `EthanolInhibition` and the two per-rate
      `ArrheniusTemperature` modifiers — the stacked config whose carbon/nitrogen
      closure is already locked in `tests/test_kinetics_arrhenius.py`. The only
      wine/beer difference is the sugar vector (1 vs 3 slots); beer's sequential uptake
      lives inside the uptake Process, so no extra Process is needed. The empty-`Medium`
      "no kinetics" baseline moved to an explicit bare `Medium` (`test_media.py`);
      `compile→simulate` now ferments end-to-end and conserves carbon for both media
      (`test_compile.py`).
- [ ] Tune functional forms + parameters against the curves.
- [ ] Unskip & pass `test_wine_24brix_ferments_to_dryness_in_10_to_14_days`.
- [ ] Unskip & pass `test_beer_1048_og_attenuates_in_5_to_7_days`.
- [ ] Unskip & pass `test_co2_integral_tracks_sugar_consumed`.
- [x] Directional check: lower T → slower (qualitative; full byproduct check is T2). →
      `test_warmer_ferments_faster_than_cooler` in `tests/test_kinetics_arrhenius.py`
      (identical uptake config, warmer run leaves less residual sugar over the span).
- [ ] Update `docs/ARCHITECTURE.md` + `DECISIONS.md`; promote tiers; commit.

Keep `pytest` / `ruff` / `mypy` green at every step. Do not weaken benchmarks.
