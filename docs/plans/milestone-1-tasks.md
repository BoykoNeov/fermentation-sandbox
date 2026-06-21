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
- [ ] `SugarUptakeToEthanolCO2` Process (Gay-Lussac yield; vector-aware for beer)
      + unit test.
- [ ] `EthanolInhibition` Process + unit test.
- [ ] **Wire parameter-tier propagation into `ProcessSet.tier_of`** — each Process
      declares the params it reads; an output's tier = lowest of (its Processes'
      tiers, those params' tiers). Closes the D-1 gap so a VALIDATED process on
      speculative params reports speculative. Add a test asserting exactly this.
- [ ] `ArrheniusTemperature` modifier hook (isothermal runs for M1) + unit test.
- [ ] Source + reconcile parameters; replace `wine_generic.yaml` placeholders and
      add `beer_generic.yaml`; re-tag tiers. (Biggest effort — handoff §2.3.)
- [ ] Tune functional forms + parameters against the curves.
- [ ] Unskip & pass `test_wine_24brix_ferments_to_dryness_in_10_to_14_days`.
- [ ] Unskip & pass `test_beer_1048_og_attenuates_in_5_to_7_days`.
- [ ] Unskip & pass `test_co2_integral_tracks_sugar_consumed`.
- [ ] Directional check: lower T → slower (qualitative; full byproduct check is T2).
- [ ] Update `docs/ARCHITECTURE.md` + `DECISIONS.md`; promote tiers; commit.

Keep `pytest` / `ruff` / `mypy` green at every step. Do not weaken benchmarks.
