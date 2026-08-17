---
name: nitrogen-dose-salt-charge
description: "D-210 - add_dap doses a SALT: both of DAP's ions are in the charge balance, and which half owns the endpoint is the opposite of what D-209 predicted"
metadata: 
  node_type: memory
  type: project
  originSessionId: 41ccdfd6-1486-4600-8cb4-43288814b4cf
  modified: 2026-08-17T13:02:01.295Z
---

**Live prohibitions — the nitrogen DOSE's charge (D-210).** Detail split out of
`.claude/memory/project-fermentation-sandbox.md`; that file's ledger points here by path. Read
it before proposing anything about `add_dap`, dosed phosphate, or the nitrogen pool's charge.
Every bullet is *what it forbids* + the record to read for *why*. **If a prohibition looks
unconvincing, go read D-210 — do not argue past it from this file.**

**`add_dap` doses a SALT and both ions SHIP. Do not re-propose either as unbuilt.**
- **`phosphate` is a state slot in BOTH registries** — `dap_phosphate_fraction` 0.74206 g
  H₃PO₄-equivalent/g DAP, VALIDATED. **Both media on purpose**: `add_dap` is medium-agnostic, so
  a wine-only scope would leave a dosed BEER booking ammonium with no anion — a strong-base
  artefact worse than the omission. **DIPROTIC on purpose** (2.15/7.20): pKa₃ 12.35 is beyond
  both the beverage and TA's 8.2 endpoint, and `AcidSpec.protons` is what TA *subtracts from*, so
  carrying it invents **0.63 g/L of TA** at 1.1 g/L DAP. **NEVER call this D-178's malt
  phosphate** — that one is present at t=0 and the anchor absorbs it; a dose lands after.
- **The dosed charge rides `nitrogen_charge_excess`, and it stores the EXCESS not the mean
  charge** — because 0.0 is then both the default and the correct undosed value, so no sentinel
  is needed and a sentinel-vs-state-float is a gate `num_jac` straddles. `dap_nitrogen_charge`
  = 1.0 VALIDATED (stoichiometry; D-209's averages are DERIVED). **One dimensionless slot, NOT a
  second N pool** — charge per mole is invariant under proportional drawdown, so no Process
  touches it and `remix_nitrogen_charge_excess` re-mixes at each dose. **The dose-time pH RISE is
  EMERGENT** (+2 cation − 0.95 anion = the protons HPO₄²⁻ takes up) — never add a proton term.
- **WHICH HALF OWNS WHAT is the opposite of D-209 §8c's sizing.** The ammonium half moves a DRY
  wine's endpoint by **exactly 0.0** (the nitrogen leaves either way) and owns the **excursion**
  (+0.235 pH at the dose vs +0.066). The **phosphate owns 100 % of the permanent** change:
  **−0.162 pH** at Palma's 1.1 g/L, **−0.043** at 0.3 g/L. The ammonium half DOES reach an
  endpoint where nitrogen is **LEFT STANDING** (+0.228 pH, late dose) — stuck ferments, sweet
  wines, short horizons. **Never quote one with/without number**: they oppose at the dose instant
  and the phosphate cancels 71 % of the bump.
- **The eight Processes that ADD to `N` keep the medium average, and that is MEASURED.** ~88 % of
  the inflow is `AminoAcidAssimilation`, an **UN-DRAW not a deamination** (composition unchanged
  ⇒ average exactly right). The rest are real deaminations whose true charge is **0 or +1 per
  channel** (excreted keto acid vs decarboxylated), bracketing the average at **≤0.033 pH**
  transiently, **~2e-5** at the endpoint. **A dose is the one channel entitled to +1** — the only
  one where the model books BOTH arriving ions.
- **Native must phosphate is out of scope on the ANCHOR, never on smallness** — sourced
  **100-500 mg/L as PO₄** (*Concepts in Wine Chemistry* Ch. 1, which also states DAP brings
  phosphate along and names H₂PO₄⁻). Absorbed at t=0; un-absorbed remainder **1.3 %** of wine's
  buffer capacity. US ceiling **960 mg/L** = 1.4-6.9× the native pool — **its unit is a named
  FORK** (mg DAP vs mg N), nothing ships from it.
- **BOTH charge writes are ATOMIC and ride D-179's gate** — found in review AFTER the first green
  suite. Ungated they were gated DIFFERENTLY: the nitrogen half checks
  `charge_balance_is_populated`, `phosphate` is a plain acid slot. On an unanchored wine the dose
  booked the ANION alone **and OPENED the gate** (acid slots are what it tests), so a *nutrient
  addition* switched the whole D-209 term on: pH 3.103 → **4.530** at the dose, **−0.647** at the
  end. **That is the Palma benchmark's shape** — its 37 stayed green because they score sugar and
  ethanol. **Never gate only one half**; the `N` jump itself is NOT gated (D-36's H₂S gate).
- **`NitrogenExceedsCationDemandError` exists because `_verb_set_ph` MASKED the diagnosis** — it
  caught bare `ValueError` and rewrote every failure as *"below the acid load's intrinsic pH"*,
  wrong in **every** dosed case probed, and its printed floor still contains the nitrogen charge.
  `set_ph` reaches D-209's guard for ~12 h after a dose (0.0485 vs 0.0345 mol⁺/L). **Never
  collapse the two branches** — a positive control pins the floor branch still fires.
- **Unbuilt and stated: diammonium SULFATE** (Handbook Vol 2 prefers it; sulfate carries ~1.96
  charges vs DAP's ~0.95, so **~2× the acidification per mole N** — no `add_das`, do not reuse
  DAP's numbers) and **the dosed `amino_acids` pool's own +1** (arginine, charge-inactive slot,
  **0.056 pH** on a maximal arm, sign genuinely unclear because autolysis REFILLS that pool).

**SUPERSEDED AT D-214 — read this before the paragraph below.** Both items are now **REFUSED**,
not parked: antiport has **zero beer-text sourcing**, and trub is a **pre-pitch** event already
inside the calibration whose post-anchor form is a **charge violation**. The paragraph below is
kept for *why they were parked in the first place*; it is no longer a reason to revisit either.
→ `.claude/memory/prohibitions/trub-settling-and-the-peptide-pair.md`

**Why D-210 did NOT take the other two items D-209 §8 named.** K⁺/H⁺ antiport and trub protein
settling both push the **same way** as D-209's term, whose high `z̄` edge already lands at pH
**4.783** against a window whose floor is **4.780** — 0.003 pH. Adding same-sign beer
acidification now walks day 7 out of Tyrell's envelope from the other side, and D-209 §7 already
located the residual as **uptake TIMING**. ~~Neither is refused — both wait on the timing beat.~~
**BOTH ARE NOW REFUSED (D-214)**: antiport on sourcing (zero beer-text hits), trub because it is a
**pre-pitch** event already inside the calibration whose post-anchor form is a charge violation.
Do not re-propose either. → `.claude/memory/prohibitions/trub-settling-and-the-peptide-pair.md`
