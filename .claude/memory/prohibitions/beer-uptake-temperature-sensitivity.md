---
name: beer-uptake-temperature-sensitivity
description: "D-217 - the corpus has NO beer sugar-uptake activation energy (searched, counted, closed), and Tyrell's trial temperature is an inference the source settles only comparatively"
metadata: 
  node_type: memory
  type: project
  originSessionId: 44f9c73d-c9c3-43fa-889b-feafa52223be
  modified: 2026-08-18T07:40:14.621Z
---

**Live prohibitions — beer's uptake temperature sensitivity and Tyrell's trial frame (D-217).**
Detail split out of `.claude/memory/project-fermentation-sandbox.md`; that file's ledger points
here by path. Read it before proposing anything about `E_a_uptake`, re-sourcing its band, or the
temperature any Tyrell comparison is scored at. Every bullet is *what it forbids* + the record to
read for *why*. **If a prohibition looks unconvincing, go read D-217 — do not argue past it from
this file.**

**MEASURED NULL and a provenance repair. No `src/` change, no parameter VALUE moved.**

- **The corpus CANNOT re-source `E_a_uptake`, and this was SEARCHED, not assumed.** 26 files (24
  texts + Peyer's thesis + Tyrell in full) × 10 patterns → **75 hits, all 75 read**; then a
  SECOND census for the other half of the question — a convertible two-temperature pair contains
  none of those patterns — 6 pattern families × the 7 beer sources → **31 hits, 28 distinct, all
  read**. **ZERO** give an activation energy, a Q10, or a convertible rate pair for wort sugar
  uptake, in beer or wine. **Do not re-run this search without a NEW text.** The band stands
  **for want of a source, not for want of looking**, and the YAML's *"pending review"* is gone.
- **Never adopt the 30,000 low edge** — still the debunked "~35 kJ/mol beer figure" residue
  (D-19, D-216 §6). Nothing found here rehabilitates it.
- **NEVER read de Andres-Toro's −97 kJ/mol as a fit, however well it lands.** −90,000 closes
  Tyrell's day 2 almost exactly and that is a **coincidence with three independent refusals**: it
  is a lumped coefficient of *their* rate law; it is the **wrong SIGN** for warmer-ferments-faster;
  and on an isothermal trial it is **perfectly degenerate with `q_sugar_max`** (reproduces D-216's
  q = 1.397 row to three decimals, overshoot on days 1 and 3 included). It is not a mechanism here.
- **The benchmark's inertness to `E_a_uptake` is EXACT, not approximate** — swept −97,000 to
  +80,000 J/mol, §2.2's days-to-FG moved **0.0000 d**. Guarded in `tests/test_organic_acids.py`
  §11. A RED there names a change to the **frame** (the benchmark leaving `T_ref`), never a rate.
- **The lever's size IS the frame's distance from `T_ref`** — the whole printed band moves
  Tyrell's day 2 by **0.0449 at 15 °C** (2.81× → 2.41×, under a fifth of the gap) and by
  **exactly 0.000000 at 20 °C**. D-216 §6's "only lever that decouples" is **conditional**, and
  had the trial run at `T_ref` the escape would not exist at all.
- **Tyrell's tube trial temperature is NOT PRINTED, and 15 °C was an inference — but it is
  CORRECT, so never "fix" it.** Tyrell ran **two** trials: §2.4.1's flasks are the ones stating
  *"isothermal at 20 °C"* and are **not** Fig. 4's source; Fig. 4 is §2.4.2's EBC tubes, whose
  only printed temperature is the **FILL** (*"cooled down to 15 °C"*). D-178 said *"pitched at"*;
  **D-207 dropped the qualifier** and nine records inherited it. §3.2 settles the direction —
  *"higher fermentation temperature … of **lab scale**"* — so the tubes ran **below 20 °C**.
  D-207 carries a ⚠ for its **GROUNDS, not its value**; nothing was retracted.
  `TYRELL_TRIAL_CELSIUS` now ships that citation.
- **The 15-20 °C ambiguity is worth 2.81× → 1.94×** of D-215's headline (~31 %), and the
  **unpaired** arm said 2.81× → 1.61× (57 %) — nearly half of the apparent frame effect is really
  the stale growth rate. **Always refit `mu_max` per frame** by D-211's own objective; the 15 °C
  refit returns 0.0340 exactly, which is the positive control.
- **The pH course CANNOT adjudicate the frame**, though an unpaired read says it can: 7/8 days
  inside at **every** temperature once `mu_max` is refit. The apparent degradation with warmth
  (5/6 → 4/6) was the stale rate. Never cite the pH as a second observable on trial temperature.
- **D-216's refusal SURVIVES the reframe** — the `q` matching Tyrell is in band and breaks the
  benchmark in **both** frames (1.397 → 2.71 d at 15 °C; 0.939 → 4.00 d at 20 °C). Only its
  MAGNITUDE was frame-dependent. **D-216 §11's anchor choice is untouched and still the owner's.**
- **A mutation arm must be IN BAND to test anything** — the store validates at load, so an
  out-of-band arm dies in pydantic and its RED says nothing. Cost two arms in D-216 and two again
  here. The §11 lever guard also fires on `E_a_growth` (its pinned 0.0449 rides the biomass);
  that is documented, not uptake-specific.
