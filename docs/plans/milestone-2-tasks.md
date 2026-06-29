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
- [x] **Carbon-accounting sub-decision SETTLED — option (b), and it is carbon-*correct*,
      not a compromise** (overrides the plan's tentative route-from-source). `Byp` already
      books "higher alcohols" (= fusels) as succinic, and ester carbon re-expresses
      ethanol+acid carbon already counted, so adding these pools to `total_carbon` would
      **double-count**; routing from sugar would force carving fusels out of `Byp`
      (re-anchoring `Y_byproduct_sugar`, risking the ABV/realised-yield guard) for zero
      benchmark gain (directional). So: pools NOT carbon-routed, `total_carbon` untouched,
      closure byte-for-byte. Boundary note recorded for any future router. Formalise in
      **D-19** at beat close.
- [ ] Source + reconcile ester/fusel rate + T-sensitivity parameters; replace
      placeholders; tag tiers honestly (`plausible`/`speculative`, directional only).
      **Hard constraint:** keep each `E_a` > `E_a_uptake` (the ordering that makes the
      run-integrated total fall with T — verified load-bearing).
- [ ] Add a multi-temperature comparison; **unskip & pass
      `test_lower_temperature_is_slower_but_cleaner`** (lower T ⇒ longer to dryness AND
      fewer esters+fusels). Confirm the §2.2 trio + carbon conservation stay green.
      **Direction already verified empirically** (scratch run): wine 14 °C→15.8 d/0.141 g/L
      vs 25 °C→5.4 d/0.186 g/L; beer 14 °C→10.8 d/0.111 g/L vs 25 °C→4.4 d/0.147 g/L —
      slower AND cleaner when colder, both media, at dryness and at run end.
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
