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
- [~] **Carbon-accounting sub-decision — DECISION: option (a), route carbon from sugar**
      (user call, 2026-06-29, overriding the advisor/author lean toward (b)). **Not yet
      implemented** — the shipped code is the **interim (b)** (esters/fusels produced-only,
      *outside* `total_carbon`). Implementing (a) is **deferred to the next session**; see
      the dedicated planned-work section below. Formalise in **D-19** once (a) lands.
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

## Planned next session — carbon-accounting option (a)

**Decision (user, 2026-06-29):** route ester/fusel carbon **from sugar** and weight
the pools in `total_carbon`, so esters/fusels become real carbon-accounted state
rather than diagnostic re-expressions. This overrides the advisor/author lean toward
(b) (the interim shipped code: pools produced-only, *outside* `total_carbon`, closure
byte-for-byte). **Not yet implemented — interim (b) stands until this lands.** The
rationale for (a): one consistent rule for every produced-only pool (`Gly`/`Byp`
already route from sugar, D-16), and machine-precision carbon closure that includes
the aroma pools rather than asserting "their carbon is booked elsewhere." Formalise
the outcome as **D-19** once (a) is in.

The hard part is the **`Byp` double-count**: `Byp` is "organic acids **+ higher
alcohols**, carbon-accounted as succinic acid" (`chemistry.M_SUCCINIC`,
`conservation.total_carbon` weights `Byp` by `carbon_mass_fraction("succinic_acid")`).
Fusels *are* higher alcohols — already inside `Byp`. Adding a separate carbon-weighted
`fusels` pool on top of the existing `Byp` weight double-books that carbon. So (a) is
not just "weight the new pools"; it forces re-anchoring what `Byp` represents.

**Implementation steps (next session):**

1. **Add representative species to `core/chemistry.py`.** Ester ⇒ ethyl acetate
   `C4H8O2` (`M ≈ 88.11`, carbon fraction 4·12.011/M); fusel ⇒ isoamyl alcohol
   `C5H12O` (`M ≈ 88.15`, carbon fraction 5·12.011/M). Add `M_ETHYL_ACETATE` /
   `M_ISOAMYL_OH` and their `MOLAR_MASS` / `CARBON_ATOMS` entries (the single source
   of truth the conservation check and kinetics both read).

2. **Routing variant — (a1) is the chosen default; (a2) is a flagged alternative.**
   - **(a1) route from sugar + carve `Byp` — THIS IS option (a) as the user chose it.**
     "Route ester/fusel carbon **from sugar** and weight the pools" *is* (a1): the
     byproduct Processes also draw the species' carbon out of `S` (touches gains `S`),
     and `Byp` is **redefined to organic-acids-only** (drop "higher alcohols") with
     `Y_byproduct_sugar` re-anchored downward so the realised-yield split still sums
     correctly. Most physically literal; highest blast radius (re-sources a Tier-2
     yield parameter, re-touches the realised-yield/ABV guard). **Implement this
     unless the user redirects.**
   - **(a2) transfer carbon from already-booked pools, no new sugar draw — a distinct
     third option to FLAG, not a flavour of (a).** Esters decrement `E` (+ a little
     `Byp`); fusels decrement `Byp` (the higher-alcohol share already there). No `S`
     touch, no `Y_byproduct_sugar` change — `total_carbon` is unchanged by
     construction and closure stays exact. **Be honest about what this is:** it does
     *not* literally route from sugar, so its outcome (exact closure, pools weighted,
     no real sugar→aroma carbon) sits functionally close to the **(b) the user
     rejected** — it is a *new* option, not (a). Its only advantage is the smaller
     blast radius (no `Y_byproduct_sugar` re-anchor, ABV guard untouched). **Surface
     it to the user for their call** (their standing "discuss disagreements"
     preference applies); do **not** silently pick it over (a1). The amino-acid-carbon
     caveat (step 5) is the strongest honest argument for considering it.

3. **Resolve the `Byp` double-count** per the chosen variant: (a1) carve higher
   alcohols out of `Byp`'s definition + re-anchor `Y_byproduct_sugar`; (a2) make the
   fusel pool a *transfer out of* `Byp` so the total is unchanged.

4. **Weight `esters`/`fusels` in `conservation.total_carbon`** (the `if "Gly"`/
   `if "Byp"` block) by their new species' carbon fractions — only meaningful once
   step 3 prevents the overlap.

5. **Caveat to flag in code + D-19:** the Ehrlich fusel pathway is *amino-acid*
   -derived, but the `N` pool carries **no carbon** in `total_carbon` (it is YAN, a
   nitrogen-only ledger). So "route fusel carbon from sugar" is already a carbon-source
   approximation — the carbon skeleton of a catabolic fusel comes from the amino acid,
   not directly from hexose. Document that the sugar (or `Byp`-transfer) source is a
   bookkeeping stand-in, not a claim about the metabolic carbon origin.

6. **Update the byproduct tests** (`tests/test_kinetics_byproducts.py`): `touches`
   assertions change (gains `S` under a1); add per-Process **carbon-balance**
   assertions (Δsugar-carbon == Δether/fusel-carbon, or pool-transfer neutrality under
   a2); keep the `falls-with-T` and isolability guards.

7. **Verify two invariants together:** (i) `total_carbon` closes to machine precision
   with byproducts on (new assertion in the conservation tests); (ii) the
   realised-yield / ABV realism guard stays green — `test_wine_abv_and_glycerol`
   asserts realised `Y_E ∈ [0.46, 0.50]` (~0.482) **and** `ABV ∈ [13.5, 15.5] %`
   (~15.0 %); neither band may move, and the glycerol band must hold too. (a1)'s
   `Y_byproduct_sugar` re-anchor is the thing most likely to nudge these; (a2) would
   leave them untouched (which is also why (a2) is closer to the rejected (b)).

**When (a) lands:** flip the carbon item above to `[x]`, rewrite the `byproducts.py`
carbon docstring (currently documents interim (b)), and record the variant chosen +
the `Byp` resolution + the amino-acid-carbon caveat in **D-19**.

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
