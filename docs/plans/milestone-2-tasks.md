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
- [ ] `FuselAlcoholsEhrlich` Process + unit test — additive, produced-only, monotone in
      T, togglable-off. Ehrlich pathway (amino-acid-derived); flag the non-monotonic-N
      simplification.
- [ ] `EsterSynthesis` Process + unit test — additive, produced-only, monotone in T
      (warmth favours esters), coupled to fermentative flux + N.
- [ ] Settle the **carbon-accounting sub-decision** (recommend route-from-source, D-16
      style → machine-precision closure); extend `total_carbon` if pools are
      carbon-routed.
- [ ] Source + reconcile ester/fusel rate + T-sensitivity parameters; replace
      placeholders; tag tiers honestly (`plausible`/`speculative`, directional only).
- [ ] Add a multi-temperature comparison; **unskip & pass
      `test_lower_temperature_is_slower_but_cleaner`** (lower T ⇒ longer to dryness AND
      fewer esters+fusels). Confirm the §2.2 trio + carbon conservation stay green.
- [ ] Record outcomes in **DECISIONS D-19**; update `milestone-2-plan.md` + ARCHITECTURE.

## Parallel (physics-free) — stochastic ensemble wrapper

- [ ] Runtime wrapper over `simulate` sampling each parameter within its provenance
      `Uncertainty` band; run ensembles; aggregate (median + spread). Lives in `runtime`,
      outside the pure core (handoff §1.6). Tier/uncertainty bands reported on outputs.
- [ ] Tests: determinism preserved for a single (unsampled) run; ensemble spread tracks
      the input uncertainty; core stays reproducible.

## Next beat — pH / acid charge-balance solver (keystone, D-18)

- [ ] Track tartaric/malic/lactic/acetic (± carbonic) as carbon-accounted state.
- [ ] Per-RHS `Σ charge = 0` root-find for `[H⁺]`; `pH` as a derived pure function
      (smooth/fast for BDF). pKa set as parameters.
- [ ] Resolve the three D-18 couplings: evolved-vs-dissolved CO₂ (carbonic); acid carbon
      vs the D-16 `Byp`=succinic sink (no double-count); pKa(T).
- [ ] Per-acid initial concentrations wired through `compile_scenario`.

## Later beats (stubs, dependency-ordered)

- [ ] **SO₂** — free/bound/molecular equilibrium, pH-driven speciation (pKa ≈ 1.81);
      acetaldehyde binding. (Needs pH.)
- [ ] **MLF (*Oenococcus oeni*)** — second-organism Process on a "pitch MLF" event;
      malic→lactic+CO₂; pH/ethanol/SO₂/T-sensitive growth. (Needs pH.)
- [ ] **Mixed cultures / Brett / sour consortium** — resource competition. (After MLF.)
- [ ] **Remaining §3.2 byproducts** — diacetyl (VDK, the lager rest), acetaldehyde
      (early transient peak), H₂S (N/S-deficiency signal).
