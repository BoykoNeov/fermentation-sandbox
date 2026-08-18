# Per-cell yeast dry mass — SETTLED at D-219, and what it forbids

**The value is 4 × 10⁻¹¹ g (40 pg) and it is NOT a literature pick.** It lives in
`fermentation/units/convert.py` as `cells_per_ml_to_pitch_gpl`, with the argument in its docstring.

## Do not re-open it as "which estimate"
Coleman, Fish & Block 2007 — the paper wine's `Y_X/N` + its N regression, `k'_d`, `mu_max` and the
`X`/`X_A` split are fitted to — states verbatim: *"Each cell count was converted to grams per liter
of cell mass, assuming that each cell weighs 4 × 10⁻¹¹ g."* He counted and weighed nothing, so
**every gram in this engine's wine biomass is a count × that constant.** It is the DEFINITION of the
unit, so a counted pitch converted at any other figure is in a unit the model's parameters do not
use. Assert it EXACTLY; an identity has no tolerance.

- **Tier is `plausible`, never `validated`** — Coleman writes *"assuming"*. Nothing in the chain is a
  weighing. Do not let a note imply otherwise.
- **Corroboration (independent, cancels his assumption):** invert to cells per g N (2.515e11 at
  330 mg N/L) and price with Roels `CH₁.₈O₀.₅N₀.₂` ⇒ **34.87 pg**, 13 % agreement. It also fixes the
  frame Coleman left open — read as WET, 40 pg makes yeast 33 % N on a dry basis (real 7-12 %).
- **Band 28-50 pg**, DERIVED live from `biomass_N_fraction`'s own 0.08-0.14, not asserted.

## Both readings the archive shipped are RETIRED — do not resurrect either
- **18 pg** was an unsourced assertion in the two wine benchmarks (now corrected to 0.0400 g/L).
  Implies a ~50 fL cell: a lab haploid, not the diploid EC1118 those files run.
- **~100 pg** was **back-computed** from `TYRELL_SCENARIO`'s `pitch_gpl = 1.0`, so it is a residual
  absorbing the cell mass **and** every per-gram rate error. Implies a ~300 fL cell. Not a cell mass.
- **Never build a band from active-dry-yeast dosing conventions.** Handbook of Enology gives
  100-200 pg (p.95) and 25-50 pg (p.334) — mutually inconsistent, and a mass of *product*
  (moisture, carrier, dead cells), not of cell. The ~10¹⁰ cells/g commercial figure is the same
  frame error and is the likeliest ORIGIN of a `pitch_gpl` of 1.0.
- Tyrell 2013 prints **no** pitch by weight (propagated, harvested at Hochkräusen, pitched by count).
  Foster 2022 measures a DCW but never pairs it with a count. Neither can settle it — checked.

## What it decided, and what is now FORBIDDEN to say
- **§2.2's 5-7 d window does NOT survive.** D-218's surviving corner was the 100 pg row; the settled
  row takes the benchmark to **3.63 d**. Verdict recorded, `q_sugar_max` NOT moved — D-216 §4's shape
  refusal stands (a Foster-satisfying rate finishes a day early against both tails and still misses
  both day 2s).
- **D-218 §4's "3.33 d, an 11 % miss … the model is not badly wrong about beer's speed" is CORRECTED.**
  At the settled conversion the shipped model takes **4.84 d against a published ≤3 (1.61×)** and
  **10.74 d against a ≤10 ceiling** at 12 °C. The other half — *the brief* is wrong — stands.
- **D-216 §7's 3.51× is CORRECTED to 5.31×** (prose only, never asserted).
- **`TYRELL_SCENARIO` carries 2.51× Tyrell's counted 0.3984 g/L.** Deliberately NOT corrected: the
  price is a `mu_max` refit (fitted at 1.0), day-2 shortfall 2.81× → **6.58×**, N drawn by 24 h
  0.360 → **0.145 (outside** Tyrell's 0.234-0.448), and beer's pH 7/8 → 6/8 days.
- **2.51× and D-215's "~2.8× too slow" are NOT two routes agreeing.** D-215's was measured at pitch
  1.0, with the excess already in. They **COMPOUND**: 6.58 ≈ 2.51 × 2.62.

## Wine call sites
Both now convert 10⁶ cells/mL → **0.0400 g/L**. **All six Varela pins held unchanged** (the repair is
nearly inert at a research pitch). **One Palma pin moved and the direction is a COST**: the
nitrogen-limited LF at 144 h goes 41.0 → 28.6 g/L against a measured ~80, so that shortfall goes
~2× → ~2.8×. Re-banded [35,48] → [24,33]. Two other Palma pins now sit near their floors (CF 128 h
vs a 125 floor; gap 1.778 vs 1.7) and were **deliberately not re-centred**.

## Still open
Nobody has paired a **count with a gravimetric dry weight in one fermentation at one timepoint** —
that is what would make this `validated`. Cramer et al. 2002 (Coleman's antecedent) not obtained;
if *they* measured it the tier rises. And the ~42 % Coleman-vs-Varela `Y_X/N` gap (D-56 finding 3) is
now a cell-number-or-cell-mass question with no data separating them — but the elemental route
forbids resolving it on mass alone (68 pg/cell would reconcile them and implies yeast at 5.2 % N).
