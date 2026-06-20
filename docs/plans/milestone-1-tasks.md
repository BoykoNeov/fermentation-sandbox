# Milestone 1 — task checklist

- [ ] Define wine `StateSchema` ({X, S(1), E, N, T, CO2}) and beer schema
      ({X, S(3), E, N, T, CO2}).
- [ ] `compile(Scenario) -> (y0, ProcessSet, params)` — the scenario→core seam,
      converting industry units to canonical at the boundary.
- [ ] Carbon + mass conservation quantity functions for the real chemistry;
      wire into runtime assertions and tests.
- [ ] `GrowthNitrogenLimited` Process + unit test.
- [ ] `SugarUptakeToEthanolCO2` Process (Gay-Lussac yield; vector-aware for beer)
      + unit test.
- [ ] `EthanolInhibition` Process + unit test.
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
