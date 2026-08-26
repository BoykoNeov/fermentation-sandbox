---
name: ensemble-anchor-reanchor
description: "D-233 - the pH anchor is re-solved per ensemble member and that is BUILT; the peptide capacity half is measured, guarded and deliberately NOT fixed"
metadata: 
  node_type: memory
  type: project
  originSessionId: e703dc5b-abaf-4a2b-a0de-748460abea2f
  modified: 2026-08-26T11:17:02.957Z
---

**Live prohibitions — the per-member pH anchor, and the peptide pair (D-233).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read it
before proposing anything about `y0` and the ensemble, `cation_charge`, `initial_ph`,
`peptide_buffer_capacity_beer`, or "the ensemble should sample X". Every bullet is *what it
forbids* + the record to read for *why*. **If a prohibition looks unconvincing, go read D-233 —
do not argue past it from this file.**

- **The anchor re-solve is BUILT — never re-propose it as unbuilt or as a known gap.**
  `simulate_ensemble` takes `y0_for_member`; `CompiledScenario.reanchor_for_member()` supplies it
  whenever the scenario gave an `initial_ph`. Members now start at the anchored pH to **2.3e-11**,
  where they used to span **5.5062-5.7778 against 5.65** (beer, worst miss 0.1438) and
  **3.4208-3.5780 against 3.50** (wine, worst 0.0792). **BOTH media, never beer-only.**
- **D-24's "`y0` is held fixed" is CORRECTED, not repealed.** The excluded axis — scenario
  INPUTS, Brix and YAN — **stands**, and this beat adds no way to vary them. The distinction is
  what a recipe **states** vs what the engine **computes from it through sampled parameters**.
  **Never cite D-233 as licence to sample a scenario input**; that still needs its own decision.
- **NEVER quote 1.287x as the band's inflation.** That number is `pKa_peptide_buffer`'s OWN
  contribution with all 82 other parameters nominal. Over the full sampled set the day-14 band
  moves **1.008x — 0.8 %**. The defect is a **PER-MEMBER trajectory error, not a band-width
  error**; worst per-member day-14 shift 0.0346 pH. **The case for the fix is the t=0 contract
  and per-member correctness, NEVER the reported spread.** I published the 29 % version first and
  retracted it in-beat [[feedback-a-one-parameter-sweep-is-not-the-band]].
- **t=0 is the airtight claim and the reason is structural.** At t=0 the state IS `y0` and the
  anchor was solved so the pH function returns `initial_ph`, so ANY nonzero t=0 spread is **100 %
  artefact** — no legitimate component to net out. At day 14 most of the spread is real physics.
  **Never argue this defect from a day-N band.**
- **The day-14 per-member shift is CONVERGED, not mesh noise** — **0.032809 pH identical to six
  decimal places at rtol 1e-6 / 1e-8 / 1e-9**, both arms on the shipped path. It HOLDS rather than
  shrinking, so "per-member trajectory error" is MEASURED. **Do not re-open it as a mesh artefact**
  [[feedback-separate-mesh-from-coupling-by-convergence]]. t=0 needs no such check — that side is a
  closed form, so 2.3e-11 is `solve_ph`'s own residual.
- **The re-anchor moves ONE slot and that is deliberate.** A full re-run of the initial builder was
  considered and **DECLINED** — it would move seeds the beat never measured (nitrogen-dependent
  yield, hop boil, every dosed inert slot). **Never "finish" it into a general y0 rebuild.**
  `cation_charge_for_ph` **already existed** (D-186, for `set_ph`) and at t=0 reduces term for term
  to the compile seam, which is what makes the nominal draw an EXACT control, not a tolerance.
- **An unsolvable member is a FAILURE, never a fallback to the nominal `y0`.** A silent fallback
  would put an unanchored member inside the reported spread with every survivorship count reading
  clean. Guarded.
- **The peptide capacity half is MEASURED and GUARDED, NOT FIXED — and that is a decision.**
  D-214's **+0.0099 is CONFIRMED at 0.0100** (re-measured, not inherited); it is **21 %** of the
  low-pKa arm's day-14 defect and **7 %** at t=0. **Never repair it by moving the BC back-solve
  into `src`**: that would make the round-trip guard — the one that FORCED the D-180 and D-181
  re-anchors by going red — compare the root-finder against itself.
- **`test_a_drawn_peptide_pka_carries_a_wort_that_is_not_peyers_1_18` PINS A DEFECT ON PURPOSE.
  A RED there means the pair was made COHERENT — do NOT revert that beat; delete the guard and say
  so in the record.**
- **BC across the pKa band is a SPAN, not a direction.** Maximal AT the nominal by construction and
  falling off on BOTH sides: **1.116059 / 1.180 / 1.145594** at 3.86 / 4.25 / 4.50. The low edge is
  worst, but **the high edge is wrong too** — never read D-214's "1.1161-1.180" one-sidedly. The
  coherent capacities (10.99 / 10.52 / 10.70 mM) match the shipped note, so that triple is current.
- **The `pKa_*` registry is sampled IFF a pH-reading Process is active** — they reach the sampled
  set through `acetaldehyde_reduction`'s `reads`, which D-160 added *because* `reads` scopes the
  sampler. That bounds the claim; the defect is absent otherwise.
- **The compile-read-AND-sampled CENSUS is NOT run and is its own beat.** D-206's
  `must_aa_fraction_methionine` and this anchor are two known members of a set nobody has
  enumerated. `reanchor_for_member` repairs **only** the anchor, so other members are still live —
  do not read this beat as having closed the class [[feedback-a-parameter-can-be-pinned-and-drawn]].
