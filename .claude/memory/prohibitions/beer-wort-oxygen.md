---
name: beer-wort-oxygen
description: "D-213 - beer's wort aeration O2 is BUILT, seeded and stripped, and DELIBERATELY INERT; the growth coupling is declined and the O2-to-acid route stays refused"
metadata: 
  node_type: memory
  type: project
  originSessionId: 319bf67b-4e83-4c11-9ffb-068bbcb68911
  modified: 2026-08-17T16:27:53.132Z
---

**Live prohibitions — beer's wort oxygen (D-213).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read
it before proposing anything about beer's dissolved oxygen, wort aeration, or an O₂ coupling in
either medium. Every bullet is *what it forbids* + the record to read for *why*. **If a
prohibition looks unconvincing, go read D-213 — do not argue past it from this file.**

**BUILT. Do not re-propose beer's dissolved oxygen as missing** — that claim was true for
~140 records and is now false.

- **`o2_wort_aeration_beer` = 6.75 mg/L, band 5.5-8.0, `plausible`; `k_o2_uptake_beer` = 0.5
  L/(g·h), band 0.2-2.0, speculative.** BOTH aeration EDGES ARE PRINTED (*Craft Beers*); the
  nominal is the **author-constructed midpoint**, so never "correct" it to a printed figure.
  **Never blend in the same passage's "4.5-6 mg O₂/L required for CCVs"** — a requirement floor
  for a different vessel, and crossing two readings of one passage is D-209's units fork.
- **The term is DELIBERATELY INERT and that is not a defect to fix.** All three O₂ consumers
  (`OxidativeAcetaldehyde`, `PhenolicBrowning`, `EllagitanninOxidation`) are **aging-gated**, so
  nothing reads the pool during fermentation. **The owner chose the beat KNOWING this**, after
  being told. What it buys: the model no longer claims a brewery ferments anaerobically from
  t=0, and a later `begin_aging` on beer cannot oxidise against a phantom 6.75 mg/L (D-212 §7's
  hazard — the sink is why seeding is safe).
- **`begin_aging` DISABLES it — the FIRST thing that verb has ever switched off. Never remove
  that.** Left enabled it competed with the aging sinks for a dosed `add_oxygen` and ate **~45 %**
  (4.38 of 8.0 mg/L). Packaging O₂ must go to the sinks whose rates are **calibrated against
  measured depletion**, not to a fermentation-phase term with an author-estimate constant. **This
  also inverts the inertness claim's basis**: inert on the DEFAULT trajectory, but an *aged* beer
  WAS measurably changed until the fix — inertness holds because the disable makes it hold. Found
  only by the FULL suite ([[feedback-full-suite-before-green]]): domain suite green, 11 failures
  outside it, 8 expectation updates + 1 `KeyError` on a reduced parameter set + this.
- **The rate law is first-order in O₂ × biomass PRESENT, never biomass FORMED.** The
  growth-coupled form was **built out and REJECTED on timing**: sources put the uptake in the
  **lag phase, BEFORE growth**, and a yield sized from yeast sterol content (~1 % dry weight,
  ~12 mol O₂/mol ergosterol ⇒ ~0.010 g O₂/g X) empties the pool at **~26 h**, contradicting
  "rapidly disappears". Realised: 6.75 → 0.91 (4 h) → 0.097 (8 h) → 0 mg/L.
- **It must NOT be a growth-Arrhenius modifier target** — it recomputes no rate, unlike
  `AceticAcidOverflow` which must be one (D-183/D-32). Adding it would CREATE the mismatch those
  modifiers exist to prevent. It rides its own tuple because its gate differs: the acids opt in
  with `initial_ph`, this is unconditional.
- **Isolability is exact at the DERIVATIVE and NOT byte-for-byte in the trajectory** — an early
  draft said byte-for-byte and was WRONG. Seeding moves `h2s` by 1.1e-9 g/L = **3.3e-4 of its
  peak, 300× the solver `rtol`**. Proved to be the adaptive MESH by **convergence, not a bound**
  (1.1e-9 → 3.4e-12 → 1.1e-13 as rtol goes 1e-6 → 1e-9 → 1e-11), because a threshold cannot
  separate mesh from coupling. **The shipped test asserts convergence.** The pool also
  undershoots to ~−1.2e-10 g/L — an `atol` artifact; pin "physically zero" against the seed.
- **The O₂→growth coupling is REFUSED at D-258 — RE-OPENED by the owner, worked, and refused on
  MEASUREMENT. Do not re-propose it, and above all do not re-propose it on "the predictions are
  reachable now": that is D-258's own §3 and it is answered.** D-213's ground (six directional
  predictions, "none reachable in the default set") really had EXPIRED — extent is scored. The
  decline survives on two better grounds. (1) **The target contradicts itself**: `mu_max` is
  fitted on a curve NORMALISED on its own peak, so a ceiling that lands Tyrell's counted fold
  puts day-1 at 0.494-0.604 against **0.235-0.448 measured on the same panel**. (2) At the only
  target left (the printed 4-5x, §below) it is a **1.076x** move needing a yield fitted 1.8-2.4x
  off its own physiological sizing — one observable, one knob, no residual (D-206 shape). What
  WOULD re-open it: a printed **residual FAN for finished beer**, or a sourced O₂-per-new-biomass
  yield. **`mu_max`'s refit is answered PER TARGET** (void at the printed one, live at the
  counted one) — do not re-measure it. → [[feedback-a-normalised-fit-couples-level-and-timing]]
- **Two D-213 objections are DEAD and must not be re-used (D-258 §5).** (a) D-213 §4's timing
  argument rejected `y·dX/dt` as an **uptake RATE LAW**; it does **not** block a cumulative
  sterol-BUDGET ceiling, which leaves the shipped rate law untouched. (b) The "a pure O₂ cap
  predicts zero growth unaerated" objection is closed BY THE CORPUS — *Craft Beers*: "in the
  absence of aeration, yeast growth is thought to be minimal or nonexistent". Also checked and
  clean: **no double-count** with `BiomassCarryingCapacity` (whose docstring names oxygen/sterol
  limitation) — it is wine-only and opt-in, "beer carrying capacity is deferred".
- **The O₂→acetate route stays REFUSED (D-212) and D-213 REINFORCES it.** Magnitude is not the
  obstacle — 8.0 mg/L = 0.25 mmol/L bounds at **≤ 30 mg/L** of acid via `y_acetaldehyde_per_o2`'s
  best case × an impossible 100 % conversion, more than the +6 to +14 needed. But that chain is
  **aging chemistry**, and in the PDH bypass **acetate is an intermediate that is CONSUMED** for
  sterol acetyl-CoA, so more O₂ plausibly means LESS acetate. **The sign is unknown.**
- **Wine is untouched** — its O₂ stays a post-ferment dose; the Process is beer-only, asserted
  through the REGISTRY, not by running a wine.
- **D-16's stale comment is CORRECTED**: beer's acid slots were described as inert *"because beer
  still has no organic-acid producer (D-16, open)"* and beer's pH as not falling. Both overtaken
  by D-181/D-183 and D-207/D-208/D-211. Third file in two beats to carry growth-conditioned or
  overtaken prose that no sweep covered.
